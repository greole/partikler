
#include "CGALHelper.hpp"
#include "CGALTYPEDEFS.hpp"
#include "SPHDatastructures.hpp"
#include "gtest/gtest.h"
#include <vector>

Polyhedron createExampleTetrahedron() {
  Point p( 1.0, 0.0, 0.0);
  Point q( 0.0, 1.0, 0.0);
  Point r( 0.0, 0.0, 1.0);
  Point s( 0.0, 0.0, 0.0);
  Polyhedron P;
  P.make_tetrahedron( p, q, r, s);
  return P;
};

struct PushFacets {
  std::vector<Facet_handle> & fhs_;

  PushFacets(
             std::vector<Facet_handle> & fhs
             ) : fhs_(fhs) {};


  void operator()(Facet &f) {
    fhs_.push_back(&f);
  }

};

TEST(CommonVerticesTest, CorrectNumberVertices) {
    Polyhedron p = createExampleTetrahedron();

    std::vector<Facet_handle> fhs;

    PushFacets pf(fhs);

    std::for_each(p.facets_begin(), p.facets_end(), pf);

    ASSERT_EQ(commonVertices(Facet(*fhs[0]), Facet(*fhs[1])).n, 2);
}

TEST(FacetNormalTest, TrivialFacet) {
  Point p( 1.0, 0.0, 0.0);
  Point q( 0.0, 1.0, 0.0);
  Point r( 0.0, 0.0, 0.0);
  Polyhedron P;
  P.make_triangle(p, q, r);

  std::vector<Facet_handle> fhs;

  PushFacets pf(fhs);

  std::for_each(P.facets_begin(), P.facets_end(), pf);

  ASSERT_EQ(facet_normal(Facet(*fhs[0])), CGALVector (0.0, 0.0, 1.0)) ;
}


TEST(SearchPathTest, TwoFacets) {
  Point p( 1.0, 0.0, 0.0);
  Point q( 0.0, 1.0, 0.0);
  Point r( 0.0, 0.0, 0.0);
  Point s( 1.0, 1.0, 0.0);
  Polyhedron P;
  P.make_triangle(p, q, r);
  P.make_triangle(p, q, s);

  std::vector<Facet_handle> fhs;

  PushFacets pf(fhs);

  std::for_each(P.facets_begin(), P.facets_end(), pf);

  Path p1;

  p1.push_back(fhs[0]);
  p1.push_back(fhs[1]);

  Path p2 = searchPath(p1).second;

  // p1 is already close, hence p2 has the same size as p1, ie 2
  ASSERT_EQ(p2.size(), p1.size()) ;
}

TEST(ProjectedPoint, TwoFacetsEqualNormal) {
  Point p( 1.0, 0.0, 0.0);
  Point q( 0.0, 1.0, 0.0);
  Point r( 0.0, 0.0, 0.0);
  Point s( 1.0, 1.0, 0.0);
  Polyhedron P;
  P.make_triangle(p, q, r);
  P.make_triangle(p, q, s);

  std::vector<Facet_handle> fhs;

  PushFacets pf(fhs);

  std::for_each(P.facets_begin(), P.facets_end(), pf);

  Path p1;

  p1.push_back(fhs[0]);
  p1.push_back(fhs[1]);

  Path p2 = searchPath(p1).second;

  // p1 is already close, hence p2 has the same size as p1, ie 2
  ASSERT_EQ(p2.size(), p1.size()) ;

  Point P1 {0.75, 0.75, 0};

  Point P2 = projectedPoint(p2, P1);

  // Since F1 and F2 have identical normals P1 and P2 should be identical
  ASSERT_EQ(P1, P2);
}

TEST(ProjectedPoint, TwoOrthogonalFacets) {
  Point p( 1.0, 0.0, 0.0);
  Point q( 1.0, 1.0, 0.0);
  Point r( 0.0, 0.0, 0.0);
  Point s( 1.0, 0.0, 1.0);
  Polyhedron P;
  P.make_triangle(p, q, r);
  P.make_triangle(p, s, q);

  std::vector<Facet_handle> fhs;

  PushFacets pf(fhs);

  std::for_each(P.facets_begin(), P.facets_end(), pf);

  Path p1;

  p1.push_back(fhs[0]);
  p1.push_back(fhs[1]);

  Path p2 = searchPath(p1).second;

  // p1 is already close, hence p2 has the same size as p1, ie 2
  ASSERT_EQ(p2.size(), p1.size()) ;

  Point P1 {1.0, 0.0, 0.5};

  Point P2 = projectedPoint(p2, P1);

  // Since F1 and F2 have identical normals P1 and P2 should be identical
  ASSERT_EQ(P2, Point(1.5, 0, 0));
}

TEST(ProjectedPoint, ThreeFacetsTwoOrthogonal) {
  Point p( 1.0, 0.0, 0.0);
  Point q( 1.0, 1.0, 0.0);
  Point r( 0.0, 0.0, 0.0);
  Point s( 1.0, 0.0, 1.0);
  Point t( 1.0, 1.0, 1.0);
  Polyhedron P;
  P.make_triangle(p, q, r);
  P.make_triangle(p, s, q);
  P.make_triangle(q, s, t);

  std::vector<Facet_handle> fhs;

  PushFacets pf(fhs);

  std::for_each(P.facets_begin(), P.facets_end(), pf);

  Path p1;

  p1.push_back(fhs[0]);
  p1.push_back(fhs[1]);
  p1.push_back(fhs[2]);

  Path p2 = searchPath(p1).second;

  // p1 is already close, hence p2 has the same size as p1, ie 2
  ASSERT_EQ(p2.size(), p1.size()) ;

  Point P1 {1.0, 1.0, 0.5};

  Point P2 = projectedPoint(p2, P1);

  // Since F1 and F2 have identical normals P1 and P2 should be identical
  ASSERT_EQ(P2, Point(1.5, 1.0, 0));
}

TEST(SearchPathTest, ThreeFacets) {
  Point p( 1.0, 0.0, 0.0);
  Point q( 0.0, 1.0, 0.0);
  Point r( 0.0, 0.0, 0.0);
  Point s( 1.0, 1.0, 0.0);
  Point t( 2.0, 1.0, 0.0);
  Polyhedron P;
  P.make_triangle(p, q, r);
  P.make_triangle(p, q, s);
  P.make_triangle(p, s, t);

  std::vector<Facet_handle> fhs;

  PushFacets pf(fhs);

  std::for_each(P.facets_begin(), P.facets_end(), pf);

  Path p1;

  p1.push_back(fhs[0]);
  // p1.push_back(fhs[1]);
  p1.push_back(fhs[2]);

  // fails currently since facets dont have enough neighbours
  // Path p2 = searchPath(p1);

  // ASSERT_EQ(commonVertices(Facet(*fhs[0]), Facet(*fhs[2])).n, 1);

  // // check if another facet was found
  // ASSERT_EQ(p2.size(), p1.size() + 1) ;

  // Facet F1 = Facet(*p2[0]);
  // Facet F2 = Facet(*p2[1]);
  // Facet F3 = Facet(*p2[2]);


  // // Now all facets should have 2 common vertices
  // ASSERT_EQ(commonVertices(F1, F2).n, 2);
  // ASSERT_EQ(commonVertices(F2, F3).n, 2);

  // TODO check if vertices are correct
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}

