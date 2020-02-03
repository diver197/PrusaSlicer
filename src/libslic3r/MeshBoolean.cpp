#include "MeshBoolean.hpp"
#include "libslic3r/TriangleMesh.hpp"
#undef PI

// Include igl first. It defines "L" macro which then clashes with our localization
#include <igl/copyleft/cgal/mesh_boolean.h>
#undef L

// CGAL headers
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Exact_integer.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

namespace Slic3r {
namespace MeshBoolean {

using MapMatrixXfUnaligned = Eigen::Map<const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor | Eigen::DontAlign>>;
using MapMatrixXiUnaligned = Eigen::Map<const Eigen::Matrix<int,   Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor | Eigen::DontAlign>>;

TriangleMesh eigen_to_triangle_mesh(const EigenMesh &emesh)
{
    auto &VC = emesh.first; auto &FC = emesh.second;
    
    Pointf3s points(size_t(VC.rows())); 
    std::vector<Vec3crd> facets(size_t(FC.rows()));
    
    for (Eigen::Index i = 0; i < VC.rows(); ++i)
        points[size_t(i)] = VC.row(i);
    
    for (Eigen::Index i = 0; i < FC.rows(); ++i)
        facets[size_t(i)] = FC.row(i);
    
    TriangleMesh out{points, facets};
    out.require_shared_vertices();
    return out;
}

EigenMesh triangle_mesh_to_eigen(const TriangleMesh &mesh)
{
    EigenMesh emesh;
    emesh.first = MapMatrixXfUnaligned(mesh.its.vertices.front().data(),
                                       Eigen::Index(mesh.its.vertices.size()),
                                       3).cast<double>();

    emesh.second = MapMatrixXiUnaligned(mesh.its.indices.front().data(),
                                        Eigen::Index(mesh.its.indices.size()),
                                        3);
    return emesh;
}

void minus(EigenMesh &A, const EigenMesh &B)
{
    auto &[VA, FA] = A;
    auto &[VB, FB] = B;
    
    Eigen::MatrixXd VC;
    Eigen::MatrixXi FC;
    igl::MeshBooleanType boolean_type(igl::MESH_BOOLEAN_TYPE_MINUS);
    igl::copyleft::cgal::mesh_boolean(VA, FA, VB, FB, boolean_type, VC, FC);
    
    VA = std::move(VC); FA = std::move(FC);
}

void minus(TriangleMesh& A, const TriangleMesh& B)
{
    EigenMesh eA = triangle_mesh_to_eigen(A);
    minus(eA, triangle_mesh_to_eigen(B));
    A = eigen_to_triangle_mesh(eA);
}

void self_union(EigenMesh &A)
{
    EigenMesh result;
    auto &[V, F] = A;
    auto &[VC, FC] = result;

    igl::MeshBooleanType boolean_type(igl::MESH_BOOLEAN_TYPE_UNION);
    igl::copyleft::cgal::mesh_boolean(V, F, Eigen::MatrixXd(), Eigen::MatrixXi(), boolean_type, VC, FC);
    
    A = std::move(result);
}

void self_union(TriangleMesh& mesh)
{
    auto eM = triangle_mesh_to_eigen(mesh);
    self_union(eM);
    mesh = eigen_to_triangle_mesh(eM);
}

namespace cgal {

namespace CGALProc    = CGAL::Polygon_mesh_processing;
namespace CGALParams  = CGAL::Polygon_mesh_processing::parameters;
using CGALException   = CGALProc::Corefinement::Self_intersection_exception;

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using _CGALMesh = CGAL::Surface_mesh<Kernel::Point_3>;

struct CGALMesh { _CGALMesh m; };

static void triangle_mesh_to_cgal(const TriangleMesh &M, _CGALMesh &out)
{
//    using Index3 = std::array<size_t, 3>;
//    auto points  = reserved_vector<_CGALMesh::Point>(M.its.vertices.size());
//    auto indices = reserved_vector<Index3>(M.its.vertices.size());
//    for (auto &v : M.its.vertices) points.emplace_back(v.x(), v.y(), v.z());
//    for (auto &_f : M.its.indices) {
//        auto f = _f.cast<size_t>();
//        indices.emplace_back(Index3{f(0), f(1), f(2)});
//    }

//    CGALProc::orient_polygon_soup(points, indices);
//    CGALProc::polygon_soup_to_polygon_mesh(points, indices, out);
    
//    // Number the faces because 'orient_to_bound_a_volume' needs a face <--> index map
//    int index = 0;
//    for (auto face : out.faces()) face = CGAL::SM_Face_index(index++);
    
//    if(CGAL::is_closed(out))
//        CGALProc::orient_to_bound_a_volume(out);
//    else
//        std::runtime_error("Mesh not watertight");
    
//     Old stuff (faster but rather unsafe)
        for (const Vec3f &v : M.its.vertices)
            out.add_vertex(_CGALMesh::Point(v.x(), v.y(), v.z()));
    
        for (const Vec3crd &face : M.its.indices) {
            auto f = face.cast<CGAL::SM_Vertex_index>();
            out.add_face(f(0), f(1), f(2));
        }    
}

static TriangleMesh cgal_to_triangle_mesh(const _CGALMesh &cgalmesh)
{
    Pointf3s points;
    std::vector<Vec3crd> facets;
    points.reserve(cgalmesh.num_vertices());
    facets.reserve(cgalmesh.num_faces());

    for (auto &vi : cgalmesh.vertices()) {
        auto &v = cgalmesh.point(vi); // Don't ask...
        points.emplace_back(v.x().interval().sup(), v.y().interval().sup(), v.z().interval().sup());
    }

    for (auto &face : cgalmesh.faces()) {
        auto    vtc = cgalmesh.vertices_around_face(cgalmesh.halfedge(face));
        int     i   = 0;
        Vec3crd trface;
        for (auto v : vtc) trface(i++) = static_cast<unsigned>(v);
        facets.emplace_back(trface);
    }
    
    TriangleMesh out{points, facets};
    out.require_shared_vertices();
    return out;
}

std::unique_ptr<CGALMesh, CGALMeshDeleter> triangle_mesh_to_cgal(const TriangleMesh &M)
{
    std::unique_ptr<CGALMesh, CGALMeshDeleter> out(new CGALMesh{});
    triangle_mesh_to_cgal(M, out->m);
    std::cout << "does self intersect out: " << CGALProc::does_self_intersect(out->m) << std::endl;
    return out;
}

TriangleMesh cgal_to_triangle_mesh(const CGALMesh &cgalmesh)
{
    return cgal_to_triangle_mesh(cgalmesh.m);
}

bool _cgal_diff(CGALMesh &A, CGALMesh &B, CGALMesh &R)
{
    const auto &p = CGALParams::throw_on_self_intersection(true);
    return CGALProc::corefine_and_compute_difference(A.m, B.m, R.m, p, p);
}

bool _cgal_union(CGALMesh &A, CGALMesh &B, CGALMesh &R)
{
    const auto &p = CGALParams::throw_on_self_intersection(true);
    return CGALProc::corefine_and_compute_union(A.m, B.m, R.m, p, p);
}

bool _cgal_intersection(CGALMesh &A, CGALMesh &B, CGALMesh &R)
{
    const auto &p = CGALParams::throw_on_self_intersection(true);
    return CGALProc::corefine_and_compute_intersection(A.m, B.m, R.m, p, p);
}

template<class Op> void _cgal_do(Op &&op, CGALMesh &A, CGALMesh &B)
{
    bool success = false;
    try {
        CGALMesh result;
        bool badA = CGALProc::does_self_intersect(A.m);
        std::cout << "does self intersect A: " << badA << std::endl;
        bool badB = CGALProc::does_self_intersect(B.m);
        std::cout << "does self intersect B: " << badB << std::endl;
        success = op(A, B, result);

        
        bool badR = CGALProc::does_self_intersect(result.m);
        std::cout << "does self intersect result: " << badR << std::endl;
//        bool bad = CGALProc::does_self_intersect(result.m);
//        std::cout << "does self intersect: " << bad << std::endl;
//        if (bad) {
//            CGALProc::remove_self_intersections(result.m);
//        }
//        std::cout << "does self intersect: " << CGALProc::does_self_intersect(result.m) << std::endl;
        A = std::move(result);
    } catch (...) {
        success = false;
    }

    if (! success)
        throw std::runtime_error("CGAL mesh boolean operation failed.");
}

void minus(CGALMesh &A, CGALMesh &B) { _cgal_do(_cgal_diff, A, B); }
void plus(CGALMesh &A, CGALMesh &B) { _cgal_do(_cgal_union, A, B); }
void intersect(CGALMesh &A, CGALMesh &B) { _cgal_do(_cgal_intersection, A, B); }

void self_union(CGALMesh &A)
{
//     _cgal_do(_cgal_union, A, A); // TODO: this is not the way
}

template<class Op> void _mesh_boolean_do(Op &&op, TriangleMesh &A, const TriangleMesh &B)
{
    CGALMesh meshA;
    CGALMesh meshB;
    triangle_mesh_to_cgal(A, meshA.m);
    triangle_mesh_to_cgal(B, meshB.m);
    
    _cgal_do(op, meshA, meshB);
    
    A = cgal_to_triangle_mesh(meshA.m);
}

void minus(TriangleMesh &A, const TriangleMesh &B)
{
    _mesh_boolean_do(_cgal_diff, A, B);
}

void plus(TriangleMesh &A, const TriangleMesh &B)
{
    _mesh_boolean_do(_cgal_union, A, B);
}

void intersect(TriangleMesh &A, const TriangleMesh &B)
{
    _mesh_boolean_do(_cgal_intersection, A, B);
}

void self_union(TriangleMesh &m)
{
//    triangle_mesh_to_cgal(m, cgalmesh.m);
    
    TriangleMeshPtrs parts = m.split();
    
    std::unique_ptr<CGALMesh, CGALMeshDeleter> cgalmeshptr;
    for (TriangleMesh *part : parts) {
        if (!part) continue;
        part->repair();
        if (!cgalmeshptr) {
            cgalmeshptr = triangle_mesh_to_cgal(*part);
            CGALProc::remove_self_intersections(cgalmeshptr->m);
            if (parts.size() == 1) plus(*cgalmeshptr, *cgalmeshptr);
        }
        else {
            CGALMesh partmesh;
            triangle_mesh_to_cgal(*part, partmesh.m);
            CGALProc::remove_self_intersections(partmesh.m);
            CGALProc::remove_self_intersections(cgalmeshptr->m);
            plus(*cgalmeshptr, partmesh);
        }
    }
    
    plus(*cgalmeshptr, *cgalmeshptr);
    CGALProc::remove_self_intersections(cgalmeshptr->m);
    
    std::cout << "does self intersect: " << CGALProc::does_self_intersect(cgalmeshptr->m) << std::endl;
    
    if (cgalmeshptr) m = cgal_to_triangle_mesh(cgalmeshptr->m);
}

bool does_self_intersect(const CGALMesh &mesh)
{
    return CGALProc::does_self_intersect(mesh.m);
}

bool does_self_intersect(const TriangleMesh &mesh)
{
    CGALMesh cgalm;
    triangle_mesh_to_cgal(mesh, cgalm.m);
    return CGALProc::does_self_intersect(cgalm.m);
}

void CGALMeshDeleter::operator()(CGALMesh *ptr)
{
    delete ptr;
}

} // namespace cgal

} // namespace MeshBoolean
} // namespace Slic3r
