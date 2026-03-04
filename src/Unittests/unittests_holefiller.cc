
#include <gtest/gtest.h>
#include <Unittests/unittests_common.hh>
#include <OpenMesh/Tools/HoleFiller/HoleFillerT.hh>
#include <OpenMesh/Tools/Utils/MeshCheckerT.hh>

namespace {

class OpenMeshHoleFiller_Triangle : public OpenMeshBase {

    protected:

        // This function is called before each test is run
        virtual void SetUp() {

            // Do some initial stuff with the member data here...
        }

        // This function is called after all tests are through
        virtual void TearDown() {

            // Do some final stuff with the member data here...
        }

    // Member already defined in OpenMeshBase
    //Mesh mesh_;
};

/*
 * ====================================================================
 * Define tests below
 * ====================================================================
 */

/*
 */
TEST_F(OpenMeshHoleFiller_Triangle,Triangle_Hole_Filling) {

  mesh_.clear();


  bool ok = OpenMesh::IO::read_mesh(mesh_, "cube_2holes.off");

  ASSERT_TRUE(ok);

  // Check setup
  EXPECT_EQ(1456u, mesh_.n_vertices() ) << "Wrong number of vertices";
  EXPECT_EQ(2864u, mesh_.n_faces() )    << "Wrong number of faces";


  // Initialize subdivider
  OpenMesh::HoleFiller::HoleFillerT<Mesh> filler(mesh_);


  // Execute the algorithm
  filler.fill_all_holes();

  if ( std::is_same<double,typename Mesh::Scalar>() ) {
      EXPECT_EQ(1504u, mesh_.n_vertices() ) << "Wrong number of vertices after smoothing?";
      EXPECT_EQ(3004u, mesh_.n_faces() )    << "Wrong number of faces after smoothing?";
  } else {
      EXPECT_EQ(1507u, mesh_.n_vertices() ) << "Wrong number of vertices after smoothing?";
      EXPECT_EQ(3010u, mesh_.n_faces() )    << "Wrong number of faces after smoothing?";
  }

}

/*
* this would crash, as of OpenMesh 11.0
*
* Basically, we want to force a configuration where the diagonal edges are non-border edges, that means there is a
* hole, but no viable triangulation exists which isn't caught internally
*
*          v0--------v1
*           | \     / |
*           |   \ /   |
*           |   / \   |
*           | /     \ |
*          v2--------v3
* and the hole filler will attempt to fill a hole in v0->v1->v2->v3
* The additional edges and vertices later are only to ensure that the diagonals are actual interior edges.
* We have to unfortunately build the halfedge links explictly, because the add_face function or a mesh loader could unpredicatably choose slightly
* different border loops which would not encounter the same problem.
*/
TEST_F(OpenMeshHoleFiller_Triangle,Triangle_Hole_Filling_ImpossibleTriangulation) {

    mesh_.clear();

    // positions don't matter its all in the topology
    const auto v0 = mesh_.add_vertex({ 0, 0, 0 });
    const auto v1 = mesh_.add_vertex({ 0, 0, 0 });
    const auto v2 = mesh_.add_vertex({ 0, 0, 0 });
    const auto v3 = mesh_.add_vertex({ 0, 0, 0 });

    const auto e0 = mesh_.new_edge(v0, v3);
    const auto e1 = mesh_.new_edge(v3, v2);
    const auto e2 = mesh_.new_edge(v2, v1);
    const auto e3 = mesh_.new_edge(v1, v0);

    // inner hole
    mesh_.set_next_halfedge_handle(e0, e1);
    mesh_.set_next_halfedge_handle(e1, e2);
    mesh_.set_next_halfedge_handle(e2, e3);
    mesh_.set_next_halfedge_handle(e3, e0);

    // outer border
    mesh_.set_next_halfedge_handle(mesh_.opposite_halfedge_handle(e0), mesh_.opposite_halfedge_handle(e3));
    mesh_.set_next_halfedge_handle(mesh_.opposite_halfedge_handle(e3), mesh_.opposite_halfedge_handle(e2));
    mesh_.set_next_halfedge_handle(mesh_.opposite_halfedge_handle(e2), mesh_.opposite_halfedge_handle(e1));
    mesh_.set_next_halfedge_handle(mesh_.opposite_halfedge_handle(e1), mesh_.opposite_halfedge_handle(e0));

    // vertex handle
    mesh_.set_halfedge_handle(v0, e0);
    mesh_.set_halfedge_handle(v1, e3);
    mesh_.set_halfedge_handle(v2, e2);
    mesh_.set_halfedge_handle(v3, e1);

    // making diagonal interior edges v0-v2 and v1-v3

    struct DiagonalInsertConfig{
        OpenMesh::VertexHandle from;
        OpenMesh::VertexHandle to;
        OpenMesh::HalfedgeHandle e_from_incoming;
        OpenMesh::HalfedgeHandle e_from_outgoing;
        OpenMesh::HalfedgeHandle e_to_incoming;
        OpenMesh::HalfedgeHandle e_to_outgoing;
    };
    std::array<DiagonalInsertConfig,2> diagonals = {{
        {
            v0,
            v2,
            mesh_.opposite_halfedge_handle(e0),
            mesh_.opposite_halfedge_handle(e3),
            mesh_.opposite_halfedge_handle(e2),
            mesh_.opposite_halfedge_handle(e1),
        },
        {
            v1,
            v3,
            mesh_.opposite_halfedge_handle(e3),
            mesh_.opposite_halfedge_handle(e2),
            mesh_.opposite_halfedge_handle(e1),
            mesh_.opposite_halfedge_handle(e0),
        }
    }};

    for (const DiagonalInsertConfig & diag : diagonals)
    {
        // not really top/bottom, but i need a variable name that is distinct enough
        const auto v_top = mesh_.add_vertex({ 0, 0, 0 });
        const auto v_bottom = mesh_.add_vertex({ 0, 0, 0 });

        const auto e_diag_center = mesh_.new_edge(diag.from, diag.to);
        const auto e_diag_center_opposite = mesh_.opposite_halfedge_handle(e_diag_center);

        const auto e_from_top = mesh_.new_edge(diag.from, v_top);
        const auto e_to_top = mesh_.new_edge(diag.to, v_top);
        const auto e_from_bottom = mesh_.new_edge(diag.from, v_bottom);
        const auto e_to_bottom = mesh_.new_edge(diag.to, v_bottom);

        const auto f_top = mesh_.new_face();
        const auto f_bottom = mesh_.new_face();

        mesh_.set_halfedge_handle(f_top, e_diag_center);
        mesh_.set_halfedge_handle(f_bottom, e_diag_center_opposite);

        mesh_.set_halfedge_handle(v_top, mesh_.opposite_halfedge_handle(e_to_top));
        mesh_.set_halfedge_handle(v_bottom, mesh_.opposite_halfedge_handle(e_from_bottom));

        // connect the two "border" quad holes up correctly

        mesh_.set_next_halfedge_handle(diag.e_from_incoming, e_from_top);
        mesh_.set_next_halfedge_handle(e_from_top, mesh_.opposite_halfedge_handle(e_to_top));
        mesh_.set_next_halfedge_handle(mesh_.opposite_halfedge_handle(e_to_top), diag.e_to_outgoing);

        mesh_.set_next_halfedge_handle(diag.e_to_incoming, e_to_bottom);
        mesh_.set_next_halfedge_handle(e_to_bottom, mesh_.opposite_halfedge_handle(e_from_bottom));
        mesh_.set_next_halfedge_handle(mesh_.opposite_halfedge_handle(e_from_bottom), diag.e_from_outgoing);

        // connect f_top
        mesh_.set_next_halfedge_handle(e_diag_center, e_to_top);
        mesh_.set_next_halfedge_handle(e_to_top, mesh_.opposite_halfedge_handle(e_from_top));
        mesh_.set_next_halfedge_handle(mesh_.opposite_halfedge_handle(e_from_top), e_diag_center);

        mesh_.set_face_handle(e_diag_center, f_top);
        mesh_.set_face_handle(e_to_top, f_top);
        mesh_.set_face_handle(mesh_.opposite_halfedge_handle(e_from_top), f_top);

        // connect f_bottom

        mesh_.set_next_halfedge_handle(e_diag_center_opposite, e_from_bottom);
        mesh_.set_next_halfedge_handle(e_from_bottom, mesh_.opposite_halfedge_handle(e_to_bottom));
        mesh_.set_next_halfedge_handle(mesh_.opposite_halfedge_handle(e_to_bottom), e_diag_center_opposite);

        mesh_.set_face_handle(e_diag_center_opposite, f_bottom);
        mesh_.set_face_handle(e_from_bottom, f_bottom);
        mesh_.set_face_handle(mesh_.opposite_halfedge_handle(e_to_bottom), f_bottom);
        ASSERT_FALSE(mesh_.is_boundary(mesh_.edge_handle(e_diag_center)));
    }

    // optional: add some additional faces around the outer border, so that we don't have a any edges that are empty on both sides. In practice, OpenMesh ensures that at least one side of the edge always has a face
    const auto v4 = mesh_.add_vertex({ 0, 0, 0 });
    const auto v5 = mesh_.add_vertex({ 0, 0, 0 });
    const auto v6 = mesh_.add_vertex({ 0, 0, 0 });
    const auto v7 = mesh_.add_vertex({ 0, 0, 0 });

    mesh_.add_face(v3, v0, v4);
    mesh_.add_face(v0, v1, v5);
    mesh_.add_face(v1, v2, v6);
    mesh_.add_face(v2, v3, v7);

    // since this is quite the construction, show here that it is actually correct, according to mesh-checker at least
    OpenMesh::Utils::MeshCheckerT<Mesh> meshChecker(mesh_);
    ASSERT_TRUE(meshChecker.check());

    // And fill hole should not crash now, please
    OpenMesh::HoleFiller::HoleFillerT<Mesh> filler(mesh_);
    filler.fill_hole(mesh_.edge_handle(e0), 1);

    EXPECT_TRUE(meshChecker.check());
}

}
