#include <catch2/catch.hpp>
#include <test_utils.hpp>

#include <libslic3r/TriangleMesh.hpp>
#include <libslic3r/MeshBoolean.hpp>

using namespace Slic3r;

TEST_CASE("Self boolean for two spheres", "[MeshBoolean]")
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

TEST_CASE("Mesh boolean for sphere and cylinder", "[MeshBoolean]")
{
    static const double R_SPHERE = 1.;
    static const double R_CYL    = .5;
    static const double H_CYL    = .5;

    TriangleMesh sphere = make_sphere(R_SPHERE);
    TriangleMesh cyl    = make_cylinder(R_CYL, H_CYL);

    double volume_shere = sphere.volume();
    double volume_cyl   = cyl.volume(); 
    
    MeshBoolean::cgal::minus(sphere, cyl);
    
    double V = sphere.volume();
    REQUIRE(V == Approx(volume_shere - volume_cyl));
}
