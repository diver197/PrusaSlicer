#ifndef libslic3r_MeshBoolean_hpp_
#define libslic3r_MeshBoolean_hpp_

#include <memory>
#include <exception>

#include <libslic3r/TriangleMesh.hpp>
#include <Eigen/Geometry>

namespace Slic3r {

namespace MeshBoolean {

using EigenMesh = std::pair<Eigen::MatrixXd, Eigen::MatrixXi>;

TriangleMesh eigen_to_triangle_mesh(const EigenMesh &emesh);
EigenMesh triangle_mesh_to_eigen_mesh(const TriangleMesh &mesh);

void minus(TriangleMesh& A, const TriangleMesh& B);
void self_union(TriangleMesh& mesh);

namespace cgal {

struct CGALMesh;

std::unique_ptr<CGALMesh> triangle_mesh_to_cgal(const TriangleMesh &M);
TriangleMesh cgal_to_triangle_mesh(const CGALMesh &cgalmesh);
    
// Do boolean mesh difference with CGAL bypassing igl.
void minus(TriangleMesh &A, const TriangleMesh &B);

// Do self union only with CGAL.
void self_union(TriangleMesh& mesh);

// does A = A - B
// CGAL takes non-const objects as arguments. I suppose it doesn't change B but
// there is no official garantee.
void minus(CGALMesh &A, CGALMesh &B);
void plus(CGALMesh &A, CGALMesh &B);
void self_union(CGALMesh &A);

bool does_self_intersect(const TriangleMesh &mesh);
bool does_self_intersect(const CGALMesh &mesh);

}

} // namespace MeshBoolean
} // namespace Slic3r
#endif // libslic3r_MeshBoolean_hpp_
