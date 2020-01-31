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
    TriangleMesh sphere = make_sphere(1.);
    TriangleMesh cyl    = make_cylinder(.5, .5);
    
    REQUIRE(!sphere.needed_repair());
    REQUIRE(cyl.is_manifold());

    double volume_shere = sphere.volume();
    double volume_cyl   = cyl.volume();
    
    SECTION("Do create hole inside the sphere") {
        TriangleMesh result = sphere;
        
        try {
            MeshBoolean::cgal::minus(result, cyl);
        } catch (...) {
            REQUIRE(false);
        }
        
        double V = result.volume();
        
        REQUIRE(V == Approx(volume_shere - volume_cyl));
        
        REQUIRE(result.is_manifold());
        CHECK(!result.needed_repair());
        
        result.WriteOBJFile("result1.obj");
    }
    
    SECTION("Remove and re-add a piece of the sphere") {
        TriangleMesh result = sphere;
        
        // Make a cylinder going all the way through the sphere at the center
        cyl = make_cylinder(.5, 2.5);
        cyl.translate(0., 0., -1.25);
        
        REQUIRE(!cyl.needed_repair());
        REQUIRE(cyl.is_manifold());
        
        TriangleMesh holed_sphere = sphere;
        TriangleMesh missing_part = sphere;
        
        try {  
            MeshBoolean::cgal::minus(holed_sphere, cyl);
            
            REQUIRE(holed_sphere.is_manifold());
            CHECK(!holed_sphere.needed_repair());
            holed_sphere.WriteOBJFile("tmp_holed_sphere.obj");
        } catch (...) {
            REQUIRE(false);
        }
        
        try {
            MeshBoolean::cgal::intersect(missing_part, cyl);
            
            REQUIRE(missing_part.is_manifold());
            CHECK(!missing_part.needed_repair());
            missing_part.WriteOBJFile("tmp_missing_part.obj");
        } catch(...) { REQUIRE(false); }
        
        try {
            result = holed_sphere;
            MeshBoolean::cgal::plus(result, missing_part);
            
            REQUIRE(result.is_manifold());
            CHECK(!result.needed_repair());
            result.WriteOBJFile("result2.obj");
        } catch (...) { REQUIRE(false); }
        
        REQUIRE(result.volume() == Approx(volume_shere));
    }
}
