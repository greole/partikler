#include "cgal/CGALHelper.hpp"
#include "cgal/CGALTYPEDEFS.hpp"
#include "core.hpp"
#include "gtest/gtest.h"
#include <vector>

TEST(FloatField, reoderVector) {

    std::vector<size_t> idx {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};

    FloatField v(10, 0.0);

    for (size_t i = 0; i < idx.size(); i++) {
        v[i] = float(i);
    }

    reorder_vector(v, idx);

    for (size_t i = 0; i < idx.size(); i++) {
        ASSERT_EQ(v[i], float(9 - i));
    }
}

TEST(FloatField, copyStdVector) {

    std::vector<float> a(2, 2.0);

    FloatField res(2, 0.0);

    res = a;

    ASSERT_EQ(res[0], 2.0);
    ASSERT_EQ(res[1], 2.0);
}

TEST(FloatField, eagerCopy) {

    FloatField a(2, 2.0);

    FloatField res(2, 0.0);

    res = a;

    ASSERT_EQ(res[0], 2.0);
    ASSERT_EQ(res[1], 2.0);
}

TEST(FloatField, sumABFloatTest) {

    FloatField a(2, 1.0);

    FloatField b(2, 2.0);

    NeighbourFieldAB n(1, {0, 1});

    FloatField res(2, 0.0);

    sum_AB_impl(res, n, A(a) + B(b));

    ASSERT_EQ(res[0], 3.0);
    ASSERT_EQ(res[1], 3.0);

    sum_AB_impl(res, n, A(a)*B(a) + A(b)*B(b));
}

TEST(FloatField, sumABVectorTest) {

    VectorField a(2, {1.0, 0, 0});

    VectorField b(2, {2.0, 0, 0});

    NeighbourFieldAB n(1, {0, 1});

    VectorField res(2, {0.0,0.0,0.0});

    sum_AB_impl(res, n, A(a) + B(b));

    ASSERT_EQ(res[0][0], 3.0);
    ASSERT_EQ(res[1][0], 3.0);

    sum_AB_impl(res, n, ab_v(a));

    ASSERT_EQ(res[0][0], 3.0); // A(a) - B(a) = 0
    ASSERT_EQ(res[1][0], 3.0);

}

TEST(FloatField, ParenthesisTest) {

    FloatField a(2, 1);
    FloatField b(2, 2);
    FloatField c(2, 3);


    FloatField res(2, 0.0);

    solve(res, (a+b)*c);

    ASSERT_EQ(res[0], 9.0);
    ASSERT_EQ(res[1], 9.0);
}

TEST(FloatField, aTest) {

    VectorField a(2, {0,0,0});
    VectorField b(2, {0,0,0});

    NeighbourFieldAB n(1, {0, 1});

    VectorField res_v(2, {0,0,0});
    FloatField res_f(2, 0);

    sum_AB_impl(res_v, n, A(a) + A(b));
    sum_AB_impl(res_f, n, (a + b)*a);
    sum_AB_impl(res_f, n, A(a)*A(b) + A(b)*A(b));

    ASSERT_EQ(res_f[0], 0.0);
    ASSERT_EQ(res_f[1], 0.0);
}

TEST(FloatField, abTest) {

    VectorField a(2, {0,0,0});
    VectorField b(2, {0,0,0});
    FloatField  c(2, 0);

    PointField p(2, {0,0,0});
    NeighbourFieldAB n(1, {0, 1});

    VectorField res_v(2, {0,0,0});
    FloatField res_f(2, 0);

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
    sum_AB_impl(res_f, n, B(a) * ab(a)) ;
    // sum_AB_impl(res_f, n, B(a) * (A(a) - B(a)) + c) ; // fails with gcc
    sum_AB_impl(res_f, n, B(a) * ab_v(a) + c) ;
    sum_AB_impl(res_f, n, ab_v(b) * ab_v(a) + ab_f(c));

    // sum_AB_impl(res_v, n, ab_v(b) * ab_p(p) );
    // sum_AB_impl(res_f, n, B(a) * (A(a) - B(a)) + c) ;

    // ASSERT_EQ(res[0][0], 0.0);
    // ASSERT_EQ(res[1][0], 0.0);
}


TEST(FloatField, abMacroTest) {

    FloatField a(2, 1.0);

    NeighbourFieldAB n(1, {0, 1});

    FloatField res(2, 0.0);

    sum_AB_impl(res, n, ab(a));

    ASSERT_EQ(res[0], 0.0);
    ASSERT_EQ(res[1], 0.0);
}

TEST(FloatField, ComplexArithmeticFields) {

    FloatField a(1, 1.0);

    FloatField b(1, 2.0);

    FloatField res(1, 0.0);

    solve(res, (a + b) / b);

    ASSERT_EQ(res[0], 1.5);
}

TEST(FloatField, MultFloatFields) {

    FloatField a(1, 1.0);

    FloatField b(1, 2.0);

    FloatField res(1, 0.0);

    solve(res, a * b);

    ASSERT_EQ(res[0], 2.0);
}

TEST(FloatField, NormFloatFields) {

    VectorField a(1, {4.0, 4.0, 2.0});

    FloatField res(1, 0.0);

    auto norm = Norm();

    solve(res, norm(a) );

    ASSERT_EQ(res[0], 6.0);
}

TEST(FloatField, PowFloatFields) {

    FloatField a(1, 2.0);

    FloatField res(1, 0.0);

    auto pow_t = Pow<float>(2.0);

    solve(res, pow_t(a) );

    ASSERT_EQ(res[0], 4.0);
}

TEST(FloatField, ManualPowFloatFields) {

    FloatField a(1, 1.0);

    FloatField res(1, 0.0);

    auto pow_t = Pow<float>(2.0);

    solve(res, a+a);

    solve(res, pow_t(res) );

    ASSERT_EQ(res[0], 4.0);
}


// TEST(FloatField, PowFloatFields2) {

//     FloatField a(1, 2.0);

//     FloatField res(1, 0.0);

//     // auto pow = Pow<float>(2.0);

//     solve(res, pow(a+a, 2.0, 1) );

//     ASSERT_EQ(res[0], 4.0);
// }

TEST(FloatField, AddFloatFields) {

    FloatField a(1, 1.0);

    FloatField b(1, 2.0);

    FloatField res(1, 0.0);

    solve(res, a + b);

    ASSERT_EQ(res[0], 3.0);
}

TEST(FloatField, AddScalarToFloatFields) {

    FloatField a(1, 1.0);

    FloatField res(1, 0.0);

    solve(res, a + 1.0);

    ASSERT_EQ(res[0], 2.0);
}

TEST(FloatField, PowFloatFieldExpr) {

    FloatField a(1, 1.0);

    FloatField b(1, 1.0);

    FloatField res(1, 0.0);

    auto pow = boost::yap::make_terminal( Pow_Wrapper<float>(2.0));

    solve(res, pow(a+b) );

    ASSERT_EQ(res[0], 4.0);
}



int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
