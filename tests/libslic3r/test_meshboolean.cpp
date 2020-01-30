#include <catch2/catch.hpp>
#include <test_utils.hpp>

#include <libslic3r/TriangleMesh.hpp>
#include <libslic3r/MeshBoolean.hpp>

using namespace Slic3r;

TEST_CASE("CGAL and TriangleMesh conversions", "[MeshBoolean]") {
    TriangleMesh sphere = make_sphere(1.);
    
    auto cgalmesh_ptr = MeshBoolean::cgal::triangle_mesh_to_cgal(sphere);
    
    REQUIRE(cgalmesh_ptr);
    
    TriangleMesh M = MeshBoolean::cgal::cgal_to_triangle_mesh(*cgalmesh_ptr);
    
    REQUIRE(M.its.vertices.size() == sphere.its.vertices.size());
    REQUIRE(M.its.indices.size() == sphere.its.indices.size());
    
    REQUIRE(M.volume() == Approx(sphere.volume()));
}

TEST_CASE("CGAL Self boolean for two spheres", "[MeshBoolean]")
{
    TriangleMesh s1 = make_sphere(1.);
    TriangleMesh s2 = make_sphere(1.);
    
    s1.translate(-0.25, 0., 0.);
    s2.translate(-0.25, 0., 0.);
    
    TriangleMesh twospheres(s1);
    twospheres.merge(s2);
    
    twospheres.require_shared_vertices();
    
    REQUIRE(MeshBoolean::cgal::does_self_intersect(twospheres));
    
    try {
        MeshBoolean::cgal::self_union(twospheres);
    } catch (...) {
        REQUIRE(false);
    }
    
    REQUIRE(! MeshBoolean::cgal::does_self_intersect(twospheres));
}

TEST_CASE("CGAL Mesh boolean for sphere and cylinder", "[MeshBoolean]")
{
    static const double R_SPHERE = 1.;
    static const double R_CYL    = .5;
    static const double H_CYL    = .5;

    TriangleMesh sphere = make_sphere(R_SPHERE);
    TriangleMesh cyl    = make_cylinder(R_CYL, H_CYL);

    double volume_shere = sphere.volume();
    double volume_cyl   = cyl.volume();
    
    try {
        MeshBoolean::cgal::minus(sphere, cyl);
    } catch (...) {
        REQUIRE(false);
    }

    double V = sphere.volume();
    REQUIRE(V == Approx(volume_shere - volume_cyl));
}
