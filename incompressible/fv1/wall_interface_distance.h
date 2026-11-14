/*
 * This file is part of the UG4 NavierStokes plugin.
 *
 * Wall and sand-air interface distance for the FV1 turbulence model.
 * Initial class template based on ShearStressFV1 and DuneNormal in shear_stress.h.
 */

#ifndef __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV1__WALL_INTERFACE_DISTANCE__
#define __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__FV1__WALL_INTERFACE_DISTANCE__

#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include "common/common.h"
#include "common/math/ugmath.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/spatial_disc/user_data/std_user_data.h"
#include "lib_disc/operator/non_linear_operator/newton_solver/newton_update_interface.h"
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_grid/tools/periodic_boundary_manager.h"
#include "lib_grid/algorithms/attachment_util.h"
#include "lib_grid/multi_grid.h"
#include "../../properties_interface.h"

#ifdef UG_PARALLEL
#include "lib_grid/parallelization/util/attachment_operations.hpp"
#include "pcl/pcl_process_communicator.h"
#endif

namespace ug{
namespace NavierStokes{

template <typename TGridFunction>
class WallInterfaceDistanceFV1
:     public StdUserData<WallInterfaceDistanceFV1<TGridFunction>, number, TGridFunction::dim>,
	  virtual public INewtonUpdate
{
	/// domain type
	typedef typename TGridFunction::domain_type domain_type;

	/// algebra type
	typedef typename TGridFunction::algebra_type algebra_type;

	/// position accessor type
	typedef typename domain_type::position_accessor_type position_accessor_type;

	/// world dimension
	static const int dim = domain_type::dim;

	/// concentration component
	static const int _C_ = domain_type::dim + 1;

	/// grid type
	typedef typename domain_type::grid_type grid_type;

	/// element type
	typedef typename TGridFunction::template dim_traits<dim>::grid_base_object elem_type;

	/// attachment accessor
	typedef PeriodicAttachmentAccessor<Vertex, ANumber> aVertexNumber;

	/// element iterator
	typedef typename TGridFunction::template dim_traits<dim>::const_iterator ElemIterator;

	/// vertex iterator
	typedef typename TGridFunction::template traits<Vertex>::const_iterator VertexIterator;

private:
	// Nodal distance cache (populated by update()).
	ANumber m_aDistance;
	aVertexNumber m_distance;

	// Mixture grid function and mesh.
	SmartPtr<TGridFunction> m_u;
	SmartPtr<ApproximationSpace<domain_type> > m_spApproxSpace;
	grid_type* m_grid;

	// Interface and update settings.
	Interface<dim>* Inter;
	bool m_distanceReady;
	
	// Reconstructed interface segments
	struct InterfaceSegment
	{
		MathVector<dim> a, b;
	};

	std::vector<InterfaceSegment> m_interfaceSegments;

	// Triangles of the piecewise-linear 3D c = interfaceValue isosurface.
	struct SurfaceTriangle { MathVector<dim> a, b, c; };
	std::vector<SurfaceTriangle> m_interfaceTriangles;

	static number dot_vec(const MathVector<dim>& a, const MathVector<dim>& b)
	{
		number out = 0.0;
		for(int i = 0; i < dim; ++i) out += a[i]*b[i];
		return out;
	}

	static MathVector<dim> subtract(const MathVector<dim>& a, const MathVector<dim>& b)
	{
		MathVector<dim> out;
		for(int i = 0; i < dim; ++i) out[i] = a[i]-b[i];
		return out;
	}

	static MathVector<dim> interpolate(const MathVector<dim>& a, const MathVector<dim>& b, number t)
	{
		MathVector<dim> out;
		for(int i = 0; i < dim; ++i) out[i] = a[i] + t*(b[i]-a[i]);
		return out;
	}

	// Ericson's closest-point triangle calculation, using UG4 MathVector.
	// Returns a nonnegative Euclidean distance. Degenerate triangles are discarded during reconstruction.
	number point_triangle_distance(const MathVector<dim>& p, const SurfaceTriangle& t) const
	{
		const MathVector<dim> ab = subtract(t.b, t.a), ac = subtract(t.c, t.a), ap = subtract(p, t.a);
		const number d1 = dot_vec(ab, ap), d2 = dot_vec(ac, ap);
		MathVector<dim> q;
		if(d1 <= 0 && d2 <= 0) q = t.a;
		else {
			const MathVector<dim> bp = subtract(p, t.b);
			const number d3 = dot_vec(ab, bp), d4 = dot_vec(ac, bp);
			if(d3 >= 0 && d4 <= d3) q = t.b;
			else {
				const number vc = d1*d4 - d3*d2;
				if(vc <= 0 && d1 >= 0 && d3 <= 0) q = interpolate(t.a, t.b, d1/(d1-d3));
				else {
					const MathVector<dim> cp = subtract(p, t.c);
					const number d5 = dot_vec(ab, cp), d6 = dot_vec(ac, cp);
					if(d6 >= 0 && d5 <= d6) q = t.c;
					else {
						const number vb = d5*d2-d1*d6;
						if(vb <= 0 && d2 >= 0 && d6 <= 0) q = interpolate(t.a, t.c, d2/(d2-d6));
						else {
							const number va = d3*d6-d5*d4;
							if(va <= 0 && d4-d3 >= 0 && d5-d6 >= 0)
								q = interpolate(t.b, t.c, (d4-d3)/((d4-d3)+(d5-d6)));
							else {
								const number denom = 1.0/(va+vb+vc);
								const number v = vb*denom, w = vc*denom;
								for(int k=0;k<dim;++k) q[k]=t.a[k]*(1.0-v-w)+t.b[k]*v+t.c[k]*w;
							}
						}
					}
				}
			}
		}
		return std::sqrt(dot_vec(subtract(p,q), subtract(p,q)));
	}

	void add_triangle(const MathVector<dim>& a, const MathVector<dim>& b, const MathVector<dim>& c)
	{
		const MathVector<dim> ab = subtract(b,a), ac = subtract(c,a);
		const number area4 = dot_vec(ab,ab)*dot_vec(ac,ac)-dot_vec(ab,ac)*dot_vec(ab,ac);
		const number scale4 = dot_vec(ab,ab)*dot_vec(ac,ac);
		if(scale4 <= 0.0 || area4 <= 1.e-24*scale4) return;
		SurfaceTriangle t; t.a=a; t.b=b; t.c=c;
		m_interfaceTriangles.push_back(t);
	}

	// Piecewise-linear interface of a P1 tetrahedron; zero values are consistently classified positive.
	void reconstruct_tetra(const MathVector<dim> x[4], const number c[4], number iso)
	{
		int pos[4], neg[4], np=0, nn=0;
		for(int i=0;i<4;++i) {
			if(c[i] >= iso) pos[np++]=i;
			else neg[nn++]=i;
		}
		if(np==0 || nn==0) return;
		const auto edge_point = [&](int a,int b) -> MathVector<dim> {
			const number t=(iso-c[a])/(c[b]-c[a]);
			return interpolate(x[a],x[b],t);
		};
		if(np==1) {
			const int a=pos[0];
			add_triangle(edge_point(a,neg[0]),edge_point(a,neg[1]),edge_point(a,neg[2]));
		} else if(nn==1) {
			const int a=neg[0];
			add_triangle(edge_point(a,pos[0]),edge_point(a,pos[1]),edge_point(a,pos[2]));
		} else {
			const MathVector<dim> a=edge_point(pos[0],neg[0]);
			const MathVector<dim> b=edge_point(pos[0],neg[1]);
			const MathVector<dim> c0=edge_point(pos[1],neg[0]);
			const MathVector<dim> d=edge_point(pos[1],neg[1]);
			add_triangle(a,b,c0);
			add_triangle(b,d,c0);
		}
	}

	
	// Acceleration only: the MPI geometry allgather above remains unchanged.
	// Leaf ranges refer to the corresponding segment/triangle vector.
	struct BVHNode
	{
		MathVector<dim> lo, hi;
		size_t begin, end;
		int left, right;
		BVHNode() : begin(0), end(0), left(-1), right(-1) {}
	};

	std::vector<BVHNode> m_bvhNodes;
	std::vector<size_t> m_bvhIndices;

	size_t num_primitives() const
	{
		return dim == 2 ? m_interfaceSegments.size() : m_interfaceTriangles.size();
	}

	void primitive_bounds(size_t i, MathVector<dim>& lo, MathVector<dim>& hi) const
	{
		if(dim == 2)
		{
			const InterfaceSegment& s = m_interfaceSegments[i];
			for(int d = 0; d < dim; ++d)
			{
				lo[d] = std::min(s.a[d], s.b[d]);
				hi[d] = std::max(s.a[d], s.b[d]);
			}
		}
		else
		{
			const SurfaceTriangle& t = m_interfaceTriangles[i];
			for(int d = 0; d < dim; ++d)
			{
				lo[d] = std::min(t.a[d], std::min(t.b[d], t.c[d]));
				hi[d] = std::max(t.a[d], std::max(t.b[d], t.c[d]));
			}
		}
	}

	number primitive_center(size_t i, int axis) const
	{
		if(dim == 2)
		{
			const InterfaceSegment& s = m_interfaceSegments[i];
			return 0.5 * (s.a[axis] + s.b[axis]);
		}
		const SurfaceTriangle& t = m_interfaceTriangles[i];
		return (t.a[axis] + t.b[axis] + t.c[axis]) / 3.0;
	}

	int build_bvh_node(size_t begin, size_t end)
	{
		const int nodeId = static_cast<int>(m_bvhNodes.size());
		m_bvhNodes.push_back(BVHNode());

		MathVector<dim> lo, hi, pLo, pHi;
		primitive_bounds(m_bvhIndices[begin], lo, hi);
		for(size_t i = begin + 1; i < end; ++i)
		{
			primitive_bounds(m_bvhIndices[i], pLo, pHi);
			for(int d = 0; d < dim; ++d)
			{
				lo[d] = std::min(lo[d], pLo[d]);
				hi[d] = std::max(hi[d], pHi[d]);
			}
		}

		m_bvhNodes[nodeId].lo = lo;
		m_bvhNodes[nodeId].hi = hi;
		m_bvhNodes[nodeId].begin = begin;
		m_bvhNodes[nodeId].end = end;

		const size_t leafSize = 8;
		if(end - begin <= leafSize) return nodeId;

		int axis = 0;
		for(int d = 1; d < dim; ++d)
			if(hi[d] - lo[d] > hi[axis] - lo[axis]) axis = d;

		const size_t mid = begin + (end - begin) / 2;
		std::nth_element(m_bvhIndices.begin() + begin,
						 m_bvhIndices.begin() + mid,
						 m_bvhIndices.begin() + end,
						 [this, axis](size_t a, size_t b) {
							 return primitive_center(a, axis) < primitive_center(b, axis);
						 });

		const int left = build_bvh_node(begin, mid);
		const int right = build_bvh_node(mid, end);
		m_bvhNodes[nodeId].left = left;
		m_bvhNodes[nodeId].right = right;
		return nodeId;
	}

	void build_bvh()
	{
		m_bvhNodes.clear();
		m_bvhIndices.clear();
		const size_t n = num_primitives();
		if(n == 0) return; // All ranks may have zero interface primitives.
		m_bvhIndices.reserve(n);
		for(size_t i = 0; i < n; ++i) m_bvhIndices.push_back(i);
		build_bvh_node(0, n); // Root is node 0.
	}

	number box_distance_squared(const MathVector<dim>& x, const BVHNode& node) const
	{
		number dist2 = 0.0;
		for(int d = 0; d < dim; ++d)
		{
			number delta = 0.0;
			if(x[d] < node.lo[d]) delta = node.lo[d] - x[d];
			else if(x[d] > node.hi[d]) delta = x[d] - node.hi[d];
			dist2 += delta * delta;
		}
		return dist2;
	}

	number nearest_bvh_node(const MathVector<dim>& x, int nodeId, number best) const
	{
		const BVHNode& node = m_bvhNodes[nodeId];
		if(box_distance_squared(x, node) >= best * best) return best;

		if(node.left < 0) // Leaf
		{
			for(size_t i = node.begin; i < node.end; ++i)
			{
				const size_t k = m_bvhIndices[i];
				const number d = (dim == 2)
					? point_segment_distance(x, m_interfaceSegments[k])
					: point_triangle_distance(x, m_interfaceTriangles[k]);
				best = std::min(best, d);
			}
			return best;
		}

		// Visit the closer child first, tightening the upper bound quickly.
		const number dl = box_distance_squared(x, m_bvhNodes[node.left]);
		const number dr = box_distance_squared(x, m_bvhNodes[node.right]);
		if(dl <= dr)
		{
			if(dl < best * best) best = nearest_bvh_node(x, node.left, best);
			if(dr < best * best) best = nearest_bvh_node(x, node.right, best);
		}
		else
		{
			if(dr < best * best) best = nearest_bvh_node(x, node.right, best);
			if(dl < best * best) best = nearest_bvh_node(x, node.left, best);
		}
		return best;
	}

	number nearest_interface_bvh(const MathVector<dim>& x, number upperBound) const
	{
		if(m_bvhNodes.empty()) return upperBound;
		return nearest_bvh_node(x, 0, upperBound);
	}

	number nearest_interface_brute_force(const MathVector<dim>& x, number upperBound) const
	{
		number best = upperBound;
		if(dim == 2)
		{
			for(size_t k = 0; k < m_interfaceSegments.size(); ++k)
				best = std::min(best, point_segment_distance(x, m_interfaceSegments[k]));
		}
		else
		{
			for(size_t k = 0; k < m_interfaceTriangles.size(); ++k)
				best = std::min(best, point_triangle_distance(x, m_interfaceTriangles[k]));
		}
		return best;
	}

	// Minimum Euclidean distance from a point to a segment
	number point_segment_distance(const MathVector<dim>& x, const InterfaceSegment& segment) const
	{
		number length2 = 0.0;
		number projection = 0.0;

		for(int d = 0; d < dim; ++d)
		{
			const number ab = segment.b[d] - segment.a[d];
			const number ax = x[d] - segment.a[d];

			length2 += ab*ab;
			projection += ax*ab;
		}

		// Handle zero-length segments
		number t = 0.0;

		if(length2 > 0.0)
			t = projection/length2;

		// Restrict projection to the segment
		if(t < 0.0) t = 0.0;
		if(t > 1.0) t = 1.0;

		// Squared distance to closest point
		number distance2 = 0.0;

		for(int d = 0; d < dim; ++d)
		{
			const number closest = segment.a[d] + t*(segment.b[d]-segment.a[d]);
			const number diff = x[d] - closest;

			distance2 += diff*diff;
		}

		return std::sqrt(distance2);
	}

public:
	/// constructor
	WallInterfaceDistanceFV1(SmartPtr<ApproximationSpace<domain_type> > approxSpace, SmartPtr<TGridFunction> spGridFct)
		: m_u(spGridFct), m_spApproxSpace(approxSpace),m_grid(NULL), Inter(NULL), m_distanceReady(false)
	{
		if(dim != 2 && dim != 3) UG_THROW("WallInterfaceDistanceFV1: Only 2D and 3D supported.");
		if(!m_u.valid()) UG_THROW("WallInterfaceDistanceFV1: Grid function is null.");
		if(!m_spApproxSpace.valid()) UG_THROW("WallInterfaceDistanceFV1: Approximation space is null.");
		if(m_u->num_fct() <= _C_) UG_THROW("WallInterfaceDistanceFV1: Concentration component is missing.");
		if(m_u->local_finite_element_id(_C_) != LFEID(LFEID::LAGRANGE, dim, 1))
			UG_THROW("WallInterfaceDistanceFV1: Concentration must use Lagrange P1 elements.");

		domain_type& domain = *m_u->domain().get();
		grid_type& grid = *domain.grid();
		m_grid = &grid;
		grid.template attach_to<Vertex>(m_aDistance);
		m_distance.access(grid, m_aDistance);
		SetAttachmentValues(m_distance, m_u->template begin<Vertex>(), m_u->template end<Vertex>(), 0.0);
	}

	virtual ~WallInterfaceDistanceFV1() {}

	void set_phase_parameters(Interface<dim>* user)
	{
		if(!user) UG_THROW("Interface pointer is null!");
		if(!user->valid())
			UG_THROW("Interface parameters have not been initialized.");

		Inter = user;
		m_distanceReady = false;
	}

	/// Update interface geometry, wall distances, and the nodal cache.
	void update()
	{
		UG_LOG("Updating Wall Interface Distance...\n");
		
		if(!Inter) UG_THROW("WallInterfaceDistanceFV1: Interface parameters have not been set.");

		
		const number interfaceValue = Inter->interface_value();

		m_distanceReady = false;
		m_interfaceSegments.clear();
		m_interfaceTriangles.clear();

		domain_type& domain = *m_u->domain().get();
		const position_accessor_type& posAcc = domain.position_accessor();

		std::vector<DoFIndex> multInd;

		for(int si = 0; si < domain.subset_handler()->num_subsets(); ++si)
		{
			ElemIterator iter = m_u->template begin<elem_type>(si);
			ElemIterator iterEnd = m_u->template end<elem_type>(si);

			for(; iter != iterEnd; ++iter)
			{
				elem_type* elem = *iter;
				const size_t numVertices = elem->num_vertices();

				if(dim == 3)
				{
					// This first stage deliberately supports only P1 tetrahedra.
					// Never silently approximate hexahedra, prisms or pyramids.
					if(elem->reference_object_id() != ROID_TETRAHEDRON)
						UG_THROW("WallInterfaceDistanceFV1: 3D currently requires tetrahedra.");
					MathVector<dim> coords3[4];
					number c3[4];
					for(size_t i=0; i<4; ++i)
					{
						Vertex* vrt=elem->vertex(i);
						coords3[i]=posAcc[vrt];
						m_u->dof_indices(vrt,_C_,multInd);
						if(multInd.size()!=1) UG_THROW("Expected one concentration DoF per vertex.");
						c3[i]=DoFRef(*m_u,multInd[0]);
					}
					reconstruct_tetra(coords3,c3,interfaceValue);
					continue;
				}

				if(numVertices != 3 && numVertices != 4)
					UG_THROW("WallInterfaceDistanceFV1: Only triangles and quadrilaterals are supported.");

				MathVector<dim> coords[4];
				number c[4];

				// Read coordinates and concentration values
				for(size_t i = 0; i < numVertices; ++i)
				{
					Vertex* vrt = elem->vertex(i);
					coords[i] = posAcc[vrt];

					m_u->dof_indices(vrt, _C_, multInd);

					if(multInd.size() != 1)
						UG_THROW("WallInterfaceDistanceFV1: Expected one concentration DoF per vertex.");

					c[i] = DoFRef(*m_u, multInd[0]);
				}
				
				// Skip elements that do not intersect the interface
				if(!Inter->cut_element(c, numVertices))
					continue;
				
				// Detect edges lying on the interface
				const number eps = 1.0e-12;
				bool alignedEdge = false;
				bool hasPositive = false;
				bool hasNegative = false;

				for(size_t i = 0; i < numVertices; ++i)
				{
					const size_t j = (i+1) % numVertices;

					const number a = c[i] - interfaceValue;
					const number b = c[j] - interfaceValue;

					if(a > eps) hasPositive = true;
					if(a < -eps) hasNegative = true;

					if(std::fabs(a) <= eps && std::fabs(b) <= eps)
						alignedEdge = true;
				}

				if(alignedEdge)
				{
					// A mixed-sign configuration may contain additional interior branches.
					if(hasPositive && hasNegative)
						UG_THROW("WallInterfaceDistanceFV1: Ambiguous contour with aligned edge.");

					if(!hasPositive && !hasNegative)
						UG_THROW("WallInterfaceDistanceFV1: Entire element lies on interface.");

					for(size_t i = 0; i < numVertices; ++i)
					{
						const size_t j = (i+1) % numVertices;

						if(std::fabs(c[i]-interfaceValue) <= eps &&
						   std::fabs(c[j]-interfaceValue) <= eps)
						{
							InterfaceSegment segment;
							segment.a = coords[i];
							segment.b = coords[j];
							m_interfaceSegments.push_back(segment);
						}
					}

					continue;
				}

				// Find intersections with c = interfaceValue
				MathVector<dim> edgePoint[4];
				std::vector<MathVector<dim> > points;

				for(size_t i = 0; i < numVertices; ++i)
				{
					const size_t j = (i+1) % numVertices;

					const number a = c[i] - interfaceValue;
					const number b = c[j] - interfaceValue;

					// Edge coincides with interface
					if(std::fabs(a) <= eps && std::fabs(b) <= eps)
						UG_THROW("WallInterfaceDistanceFV1: Interface coincides with mesh edge.");

					// Find crossing on element edge
					if((a < -eps && b > eps) || (a > eps && b < -eps) ||
					   std::fabs(a) <= eps || std::fabs(b) <= eps)
					{
						const number t = std::fabs(a) <= eps ? 0.0 :
										 (std::fabs(b) <= eps ? 1.0 : -a/(b-a));

						MathVector<dim> x;

						for(int d = 0; d < dim; ++d)
							x[d] = coords[i][d] + t*(coords[j][d]-coords[i][d]);

						edgePoint[i] = x;

						bool duplicate = false;

						for(size_t k = 0; k < points.size(); ++k)
						{
							const number dx = x[0] - points[k][0];
							const number dy = x[1] - points[k][1];

							if(dx*dx + dy*dy <= eps*eps)
							{
								duplicate = true;
								break;
							}
						}

						if(!duplicate)
							points.push_back(x);
					}
				}

				// Ordinary interface: triangle or quadrilateral
				if(points.size() == 2)
				{
					InterfaceSegment segment;
					segment.a = points[0];
					segment.b = points[1];

					m_interfaceSegments.push_back(segment);
				}
				else if(points.size() == 1)
				{
					// Isolated point on the interface
					InterfaceSegment segment;
					segment.a = points[0];
					segment.b = points[0];

					m_interfaceSegments.push_back(segment);
				}
				// Ambiguous quadrilateral: four edge crossings
				else if(numVertices == 4 && points.size() == 4)
				{
					const number f0 = c[0] - interfaceValue;
					const number f1 = c[1] - interfaceValue;
					const number f2 = c[2] - interfaceValue;
					const number f3 = c[3] - interfaceValue;

					const number D = f0*f2 - f1*f3;

					if(std::fabs(D) > eps)
					{
						InterfaceSegment s1, s2;

						if(D > 0.0)
						{
							s1.a = edgePoint[0];
							s1.b = edgePoint[1];

							s2.a = edgePoint[2];
							s2.b = edgePoint[3];
						}
						else
						{
							s1.a = edgePoint[0];
							s1.b = edgePoint[3];

							s2.a = edgePoint[1];
							s2.b = edgePoint[2];
						}

						m_interfaceSegments.push_back(s1);
						m_interfaceSegments.push_back(s2);
					}
					else
					{
						// Bilinear saddle configuration
						const number b = f1 - f0;
						const number c = f3 - f0;
						const number d = f0 - f1 + f2 - f3;

						if(std::fabs(d) <= eps)
							UG_THROW("WallInterfaceDistanceFV1: Degenerate bilinear interface.");

						const number xi = -c/d;
						const number eta = -b/d;

						MathVector<dim> saddle;

						for(int k = 0; k < dim; ++k)
							saddle[k] = (1.0-xi)*(1.0-eta)*coords[0][k] +
										xi*(1.0-eta)*coords[1][k] +
										xi*eta*coords[2][k] +
										(1.0-xi)*eta*coords[3][k];

						for(size_t i = 0; i < 4; ++i)
						{
							InterfaceSegment segment;
							segment.a = edgePoint[i];
							segment.b = saddle;

							m_interfaceSegments.push_back(segment);
						}
					}
				}
				else if(points.size() > 2)
					UG_THROW("WallInterfaceDistanceFV1: Degenerate interface configuration.");

			}
		}

		UG_LOG("WallInterfaceDistanceFV1: Reconstructed " << m_interfaceSegments.size() << " segments and " << m_interfaceTriangles.size() << " triangles (local).\n");

#ifdef UG_PARALLEL
		// Correctness-first MPI reference: exchange all reconstructed interface primitives.
		// All ranks must call update() together, even ranks without a cut element.
		if(pcl::NumProcs() > 1)
		{
			const int verticesPerPrimitive = (dim == 2 ? 2 : 3);
			const int stride = verticesPerPrimitive * dim;
			std::vector<number> localCoords, globalCoords;
			if(dim == 2)
			{
				localCoords.reserve(m_interfaceSegments.size() * stride);
				for(size_t k = 0; k < m_interfaceSegments.size(); ++k)
				{
					for(int d = 0; d < dim; ++d) localCoords.push_back(m_interfaceSegments[k].a[d]);
					for(int d = 0; d < dim; ++d) localCoords.push_back(m_interfaceSegments[k].b[d]);
				}
			}
			else
			{
				localCoords.reserve(m_interfaceTriangles.size() * stride);
				for(size_t k = 0; k < m_interfaceTriangles.size(); ++k)
				{
					for(int d = 0; d < dim; ++d) localCoords.push_back(m_interfaceTriangles[k].a[d]);
					for(int d = 0; d < dim; ++d) localCoords.push_back(m_interfaceTriangles[k].b[d]);
					for(int d = 0; d < dim; ++d) localCoords.push_back(m_interfaceTriangles[k].c[d]);
				}
			}

			pcl::ProcessCommunicator communicator;
			communicator.allgatherv(globalCoords, localCoords);
			if(globalCoords.size() % stride != 0)
				UG_THROW("WallInterfaceDistanceFV1: Invalid gathered interface coordinates.");

			if(dim == 2)
			{
				m_interfaceSegments.clear();
				for(size_t i = 0; i < globalCoords.size(); i += stride)
				{
					InterfaceSegment seg;
					for(int d = 0; d < dim; ++d)
					{
						seg.a[d] = globalCoords[i+d];
						seg.b[d] = globalCoords[i+dim+d];
					}
					m_interfaceSegments.push_back(seg);
				}
			}
			else
			{
				m_interfaceTriangles.clear();
				for(size_t i = 0; i < globalCoords.size(); i += stride)
				{
					SurfaceTriangle tri;
					for(int d = 0; d < dim; ++d)
					{
						tri.a[d] = globalCoords[i+d];
						tri.b[d] = globalCoords[i+dim+d];
						tri.c[d] = globalCoords[i+2*dim+d];
					}
					m_interfaceTriangles.push_back(tri);
				}
			}
			UG_LOG("WallInterfaceDistanceFV1: Gathered " << globalCoords.size()/stride
				   << " interface primitives on this rank.\n");
		}
#endif

		// Build one BVH over the gathered interface geometry on each MPI rank.
		// For the first comparison on a small mesh, change this flag to true.
		const bool verifyBVH = false;
		build_bvh();
		UG_LOG("WallInterfaceDistanceFV1: BVH nodes = " << m_bvhNodes.size()
			   << ", primitives = " << num_primitives() << ".\n");

		VertexIterator iter = m_u->template begin<Vertex>();
		VertexIterator iterEnd = m_u->template end<Vertex>();
		number maxVerificationError = 0.0;
		for(; iter != iterEnd; ++iter)
		{
			Vertex* vrt = *iter;
			const MathVector<dim>& x = posAcc[vrt];
			const number wallDistance = std::max((number)0.0, x[dim-1]);
			// The analytical bottom-wall distance is a rigorous upper bound
			// for the combined wall/interface distance.
			const number dist = nearest_interface_bvh(x, wallDistance);
			m_distance[vrt] = dist;

			if(verifyBVH)
			{
				const number reference = nearest_interface_brute_force(x, wallDistance);
				const number err = std::fabs(reference - dist);
				maxVerificationError = std::max(maxVerificationError, err);
				const number tol = 1.e-10 * std::max((number)1.0, wallDistance);
				if(err > tol)
					UG_THROW("WallInterfaceDistanceFV1: BVH / brute-force mismatch.");
			}
		}
		if(verifyBVH)
			UG_LOG("WallInterfaceDistanceFV1: BVH verification max abs error = "
				   << maxVerificationError << ".\n");

		UG_LOG("WallInterfaceDistanceFV1: Interface distance field updated.\n");

		m_distanceReady = true;
	}

	template <int refDim>
	inline void evaluate(number vValue[], const MathVector<dim> vGlobIP[], number time, int si, GridObject* elem, const MathVector<dim> vCornerCoords[], const MathVector<refDim> vLocIP[], const size_t nip, LocalVector* u, const MathMatrix<refDim, dim>* vJT = NULL) const
	{
		if(!m_distanceReady) UG_THROW("WallInterfaceDistanceFV1: Call update() before evaluating distance.");
		UG_ASSERT(dynamic_cast<elem_type*>(elem) != NULL, "Unsupported element type");
		elem_type* element = static_cast<elem_type*>(elem);
		const ReferenceObjectID roid = elem->reference_object_id();
		const size_t numVertices = element->num_vertices();

		const LocalShapeFunctionSet<refDim>& rTrialSpace = LocalFiniteElementProvider::get<refDim>(roid, LFEID(LFEID::LAGRANGE, refDim, 1));
		std::vector<number> shapes;
		for(size_t ip = 0; ip < nip; ++ip)
		{
			rTrialSpace.shapes(shapes, vLocIP[ip]);
			number value = 0.0;
			for(size_t sh = 0; sh < numVertices; ++sh)
				value += m_distance[element->vertex(sh)] * shapes[sh];
			vValue[ip] = value;
		}
	}

	virtual void operator()(number& value, const MathVector<dim>& globIP, number time, int si) const
	{
		UG_THROW("WallInterfaceDistanceFV1: Need element for evaluation.");
	}

	virtual void operator()(number vValue[], const MathVector<dim> vGlobIP[], number time, int si, const size_t nip) const
	{
		UG_THROW("WallInterfaceDistanceFV1: Need element for evaluation.");
	}

	virtual void compute(LocalVector* u, GridObject* elem, const MathVector<dim> vCornerCoords[], bool bDeriv = false)
	{
		const int si = this->subset();
		for(size_t s = 0; s < this->num_series(); ++s)
			evaluate<dim>(this->values(s), this->ips(s), this->time(s), si, elem, vCornerCoords, this->template local_ips<dim>(s), this->num_ip(s), u);
	}

	virtual void compute(LocalVectorTimeSeries* u, GridObject* elem, const MathVector<dim> vCornerCoords[], bool bDeriv = false)
	{
		const int si = this->subset();
		for(size_t s = 0; s < this->num_series(); ++s)
			evaluate<dim>(this->values(s), this->ips(s), this->time(s), si, elem, vCornerCoords, this->template local_ips<dim>(s), this->num_ip(s), &(u->solution(this->time_point(s))));
	}

	virtual bool continuous() const {return true;}
	virtual bool requires_grid_fct() const {return true;}
};

} // namespace NavierStokes
} // namespace ug

#endif
