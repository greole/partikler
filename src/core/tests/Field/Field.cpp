#include "cgal/CGALHelper.hpp"
#include "cgal/CGALTYPEDEFS.hpp"
#include "core.hpp"
#include "gtest/gtest.h"
#include <vector>

TEST(ScalarField, reoderVector) {

    std::vector<size_t> idx {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};

    ScalarField v(10, 0.0);

    for (size_t i = 0; i < idx.size(); i++) {
        v[i] = float(i);
    }

    reorder_vector(v, idx);

    for (size_t i = 0; i < idx.size(); i++) {
        ASSERT_EQ(v[i], float(9 - i));
    }
}

TEST(ScalarField, copyStdVector) {

    std::vector<float> a(2, 2.0);

    ScalarField res(2, 0.0);

    res = a;

    ASSERT_EQ(res[0], 2.0);
    ASSERT_EQ(res[1], 2.0);
}

TEST(ScalarField, eagerCopy) {

    ScalarField a(2, 2.0);

    ScalarField res(2, 0.0);

    res = a;

    ASSERT_EQ(res[0], 2.0);
    ASSERT_EQ(res[1], 2.0);
}

TEST(ScalarField, sumABFloatTest) {

    ScalarField a(2, 1.0);

    ScalarField b(2, 2.0);

    NeighbourFieldAB n(1, {0, 1});

    ScalarField res(2, 0.0);

    sum_AB_impl(res, n, A(a) + B(b));

    ASSERT_EQ(res[0], 3.0);
    ASSERT_EQ(res[1], 3.0);

    sum_AB_impl(res, n, A(a) * B(a) + A(b) * B(b));
}

TEST(ScalarField, sumABVectorTest) {

    VectorField a(2, {1.0, 0, 0});

    VectorField b(2, {2.0, 0, 0});

    NeighbourFieldAB n(1, {0, 1});

    VectorField res(2, {0.0, 0.0, 0.0});

    sum_AB_impl(res, n, A(a) + B(b));

    ASSERT_EQ(res[0][0], 3.0);
    ASSERT_EQ(res[1][0], 3.0);

    sum_AB_impl(res, n, ab(a));

    ASSERT_EQ(res[0][0], 3.0); // A(a) - B(a) = 0
    ASSERT_EQ(res[1][0], 3.0);
}

TEST(ScalarField, ParenthesisTest) {

    ScalarField a(2, 1);
    ScalarField b(2, 2);
    ScalarField c(2, 3);

    ScalarField res(2, 0.0);

    solve_impl(res, (a + b) * c);

    ASSERT_EQ(res[0], 9.0);
    ASSERT_EQ(res[1], 9.0);
}

TEST(ScalarField, aTest) {

    VectorField a(2, {0, 0, 0});
    VectorField b(2, {0, 0, 0});

    NeighbourFieldAB n(1, {0, 1});

    VectorField res_v(2, {0, 0, 0});
    ScalarField res_f(2, 0);

    sum_AB_impl(res_v, n, A(a) + A(b));
    sum_AB_impl(res_f, n, (a + b) * a);
    sum_AB_impl(res_f, n, A(a) * A(b) + A(b) * A(b));

    ASSERT_EQ(res_f[0], 0.0);
    ASSERT_EQ(res_f[1], 0.0);
}

TEST(ScalarField, abTest) {

    VectorField a(2, {0, 0, 0});
    VectorField b(2, {0, 0, 0});
    ScalarField c(2, 0);

    PointField p(2, {0, 0, 0});
    NeighbourFieldAB n(1, {0, 1});

    VectorField res_v(2, {0, 0, 0});
    ScalarField res_f(2, 0);

    sum_AB_impl(res_v, n, A(a) + A(b));
    sum_AB_impl(res_v, n, A(a) + (A(a) + A(b)));
    // sum_AB_impl(res, n, (A(a) + A(b)) + A(a)); // doesn't work
    // sum_AB_impl(res, n, ab(A(a)) + A(a));      // doesn't work
    // sum_AB_impl(res, n, A(a) + ab(A(a)));      // doesn't work
    // sum_AB_impl(res, n, (A(a) + A(b))*A(b));   // doesn't work

    // sum_AB_impl(res, n, ab(a) + ab(b));

    sum_AB_impl(res_f, n, A(a) * (A(a) + A(b)));
    sum_AB_impl(res_f, n, A(a) * (A(a) + B(b)));
    sum_AB_impl(res_f, n, A(a) * (A(a) + B(a)));
    sum_AB_impl(res_f, n, B(a) * (A(a) + B(a)));
    sum_AB_impl(res_f, n, B(a) * ab(a));
    // sum_AB_impl(res_f, n, B(a) * (A(a) - B(a)) + c) ; // fails with gcc
    sum_AB_impl(res_f, n, B(a) * ab_v(a) + c);
    sum_AB_impl(res_f, n, ab_v(b) * ab_v(a) + ab_f(c));

    // sum_AB_impl(res_v, n, ab_v(b) * ab_p(p) );
    // sum_AB_impl(res_f, n, B(a) * (A(a) - B(a)) + c) ;

    // ASSERT_EQ(res[0][0], 0.0);
    // ASSERT_EQ(res[1][0], 0.0);
}

TEST(ScalarField, abMacroTest) {

    ScalarField a(2, 1.0);

    NeighbourFieldAB n(1, {0, 1});

    ScalarField res(2, 0.0);

    sum_AB_impl(res, n, ab(a));

    ASSERT_EQ(res[0], 0.0);
    ASSERT_EQ(res[1], 0.0);
}

TEST(ScalarField, ComplexArithmeticFields) {

    ScalarField a(1, 1.0);

    ScalarField b(1, 2.0);

    ScalarField res(1, 0.0);

    solve_impl(res, (a + b) / b);

    ASSERT_EQ(res[0], 1.5);
}

TEST(ScalarField, MultScalarFields) {

    ScalarField a(1, 1.0);

    ScalarField b(1, 2.0);

    ScalarField res(1, 0.0);

    solve_impl(res, a * b);

    ASSERT_EQ(res[0], 2.0);
}

TEST(ScalarField, NormScalarFields) {

    VectorField a(1, {4.0, 4.0, 2.0});

    ScalarField res(1, 0.0);

    auto norm = Norm();

    solve_impl(res, norm(a));

    ASSERT_EQ(res[0], 6.0);
}

TEST(ScalarField, PowScalarFields) {

    ScalarField a(1, 2.0);

    ScalarField res(1, 0.0);

    auto pow_t = Pow<float>(2.0);

    solve_impl(res, pow_t(a));

    ASSERT_EQ(res[0], 4.0);
}

TEST(ScalarField, ManualPowScalarFields) {

    ScalarField a(1, 1.0);

    ScalarField res(1, 0.0);

    auto pow_t = Pow<float>(2.0);

    solve_impl(res, a + a);

    solve_impl(res, pow_t(res));

    ASSERT_EQ(res[0], 4.0);
}

// TEST(ScalarField, PowScalarFields2) {

//     ScalarField a(1, 2.0);

//     ScalarField res(1, 0.0);

//     // auto pow = Pow<float>(2.0);

//     solve(res, pow(a+a, 2.0, 1) );

//     ASSERT_EQ(res[0], 4.0);
// }

TEST(ScalarField, AddScalarFields) {

    ScalarField a(1, 1.0);

    ScalarField b(1, 2.0);

    ScalarField res(1, 0.0);

    solve_impl(res, a + b);

    ASSERT_EQ(res[0], 3.0);
}

TEST(ScalarField, AddScalarToScalarFields) {

    ScalarField a(1, 1.0);

    ScalarField res(1, 0.0);

    solve_impl(res, a + 1.0);

    ASSERT_EQ(res[0], 2.0);
}

TEST(ScalarField, PowScalarFieldExpr) {

    ScalarField a(1, 1.0);

    ScalarField b(1, 1.0);

    ScalarField res(1, 0.0);

    auto pow = boost::yap::make_terminal(Pow_Wrapper<float>(2.0));

    solve_impl(res, pow(a + b));

    ASSERT_EQ(res[0], 4.0);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
