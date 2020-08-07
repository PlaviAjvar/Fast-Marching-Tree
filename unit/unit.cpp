#include "gtest/gtest.h"
#include "geometry.h"
#include <vector>

/*************************
point::operator==
*************************/

TEST (equals, equals) {
    const double eps = 8e-7;
    EXPECT_TRUE(point <double>(std::vector <double>{0,0}) == point2d <double>(0,0));
    EXPECT_TRUE(point2d <double>(0,0) == point2d <double>(eps, eps));
    EXPECT_TRUE(point3d <double>(2,3,4) == point3d <double>(2+eps, 3-eps, 4+eps));

    EXPECT_FALSE(point3d <double>(2,3,4) == point3d <double>(2+2*eps, 3-2*eps, 4+2*eps));
    EXPECT_FALSE(point3d <double>(2,3,4) == point2d <double>(2,3));
    EXPECT_FALSE(point3d <double>(2,2,2) == point2d <double>(2,2));
}

/*************************
point::normalize
*************************/

// test if normalize works for nonzero vectors
TEST (normalize, nonzero) {
    point <double> P(std::vector <double>{0.1});
    point <double> Q(std::vector <double>{-1, 2});
    point <double> R(std::vector <double>{-1, 1, -1});
    point <double> S(std::vector <double>{1, -1, 1});

    point <double> Pr(std::vector <double>{1});
    point2d <double> Qr(-1.0 / sqrt(5), 2.0 / sqrt(5));
    point <double> Rr(std::vector <double>{-0.57735026, 0.57735026, -0.57735026});
    point <double> Sr(std::vector <double>{0.57735026, -0.57735026, 0.57735026});

    // test equality
    EXPECT_TRUE(P.normalize() == Pr);
    EXPECT_TRUE(Q.normalize() == Qr);
    EXPECT_TRUE(R.normalize() == Rr);
    EXPECT_TRUE(S.normalize() == Sr);

    // test inequality
    EXPECT_FALSE(P.normalize() == Qr);
    EXPECT_FALSE(P.normalize() == point <double>(std::vector <double>{-1}));
    EXPECT_FALSE(R.normalize() == Sr);
    EXPECT_FALSE(S.normalize() == Rr);
}

// test normalize for zero vector
TEST (normalize, zero) {
    point <double> P(std::vector <double>{0});
    point <double> Q(std::vector <double>{0, -0});
    point <double> R(std::vector <double>{1e-7, 1e-8, 1e-9});

    // check magnitude
    const double eps = 1e-6;
    EXPECT_TRUE(magnitude(P) < eps);
    EXPECT_TRUE(magnitude(Q) < eps);
    EXPECT_TRUE(magnitude(R) < eps);

    // normalize to itself
    EXPECT_TRUE(P.normalize() == P);
    EXPECT_TRUE(Q.normalize() == Q);
    EXPECT_TRUE(R.normalize() == R);
}

/*************************
point::operator*
*************************/

// test dot product on coplanar cases (standard)
TEST (dotproduct, coplanar) {
    point2d <double> u(1, 0);
    point2d <double> v(sqrt(2), sqrt(2));
    point3d <double> r(1, 2, 3);
    point3d <double> s(4, 1, 3);
    point3d <double> a(-1, 2, 2);
    point3d <double> b(-3, -2, 1);

    const double epsilon = 1e-6;

    EXPECT_DOUBLE_EQ(u * v, sqrt(2));
    EXPECT_DOUBLE_EQ(r * s, 15);
    EXPECT_DOUBLE_EQ(a * b, 1);

    // cos(phi) < 1
    EXPECT_TRUE(std::abs(u * v) < magnitude(u) * magnitude(v));
    EXPECT_TRUE(std::abs(r * s) < magnitude(r) * magnitude(s));
    EXPECT_TRUE(std::abs(a * b) < magnitude(a) * magnitude(b));
}

// test dot product on colinear cases
TEST (dotproduct, colinear) {
    point2d <double> u(1, 0);
    point2d <double> v(5, 0);
    point3d <double> r(1, 2, 3);
    point3d <double> s(-2, -4, -6);
    point3d <double> a(0, 0, 0);
    point3d <double> b(1, -2, -2);

    const double epsilon = 1e-6;
    EXPECT_DOUBLE_EQ(u * v, 5);
    EXPECT_DOUBLE_EQ(r * s, -28);
    EXPECT_DOUBLE_EQ(a * b, 0);

    // cos phi = 1
    EXPECT_NEAR(std::abs(u * v), magnitude(u) * magnitude(v), epsilon);
    EXPECT_NEAR(std::abs(r * s), magnitude(r) * magnitude(s), epsilon);
}

// test dot product on perpendicular cases
TEST (dotproduct, perpendicular) {
    point2d <double> u(1, 0);
    point2d <double> v(0, 2);
    point3d <double> r(1, 2, 0);
    point3d <double> s(-2, 1, 0);
    point3d <double> a(0, 0, 0);
    point3d <double> b(1, -2, -2);

    const double epsilon = 1e-6;
    EXPECT_NEAR(u * v, 0, epsilon);
    EXPECT_NEAR(r * s, 0, epsilon);
    EXPECT_NEAR(a * b, 0, epsilon);
}

// size mismatch
TEST (dotproduct, size_mismatch) {
    point3d <double> u(1, 1, 1);
    point2d <double> v(0, 0);
    point2d <double> r(1, 1);
    point <double> s(std::vector <double>{0});

    EXPECT_THROW(u * v, std::length_error);
    EXPECT_THROW(r * s, std::length_error);
}

/*************************
point3d::operator^
*************************/

// test cross product on coplanar cases (standard)
TEST (crossproduct, coplanar) {
    point3d <double> u(1, 0, 0);
    point3d <double> v(sqrt(2), sqrt(2), 0);
    point3d <double> r(1, 2, 3);
    point3d <double> s(4, 1, 3);
    point3d <double> a(-1, 2, 2);
    point3d <double> b(-3, -2, 1);

    point3d <double> uv(0, 0, sqrt(2));
    point3d <double> rs(3, 9, -7);
    point3d <double> ab(6, -5, 8);
    EXPECT_TRUE((u ^ v) == uv);
    EXPECT_TRUE((r ^ s) == rs);
    EXPECT_TRUE((a ^ b) == ab);

    // sin(phi) < 1
    EXPECT_TRUE(magnitude(u ^ v) < magnitude(u) * magnitude(v));
    EXPECT_TRUE(magnitude(r ^ s) < magnitude(r) * magnitude(s));
    EXPECT_TRUE(magnitude(a ^ b) < magnitude(a) * magnitude(b));
}

// test cross product on colinear cases
TEST (crossproduct, colinear) {
    point3d <double> u(1, 0, 1);
    point3d <double> v(5, 0, 5);
    point3d <double> r(1, 2, 3);
    point3d <double> s(-2, -4, -6);
    point3d <double> a(0, 0, 0);
    point3d <double> b(1, -2, -2);

    point3d <double> z(0, 0, 0);
    EXPECT_TRUE((u ^ v) == z);
    EXPECT_TRUE((r ^ s) == z);
    EXPECT_TRUE((a ^ b) == z);
}

// test cross product on perpendicular cases
TEST (crossproduct, perpendicular) {
    point3d <double> u(1, 0, 0);
    point3d <double> v(0, 2, 0);
    point3d <double> r(1, 2, 0);
    point3d <double> s(-2, 1, 0);
    point3d <double> a(0, 0, 0);
    point3d <double> b(1, -2, -2);

    point3d <double> uv(0, 0, 2);
    point3d <double> rs(0, 0, 5);
    point3d <double> ab(0, 0, 0);

    EXPECT_TRUE((u ^ v) == uv);
    EXPECT_TRUE((r ^ s) == rs);
    EXPECT_TRUE((a ^ b) == ab);

    // sin(phi) = 1
    const double epsilon = 1e-6;
    EXPECT_NEAR(magnitude(u ^ v), magnitude(u) * magnitude(v), epsilon);
    EXPECT_NEAR(magnitude(r ^ s), magnitude(r) * magnitude(s), epsilon);
    EXPECT_NEAR(magnitude(a ^ b), magnitude(a) * magnitude(b), epsilon);
}

/*************************
weighed_euclidean
*************************/

// test normal cases (3D)
TEST (weighedEuclidean, 3D) {
    std::vector <double> link_length{5, 5, 5};
    point3d <double> A(0, 0, 0);
    point3d <double> B(0, 0, 5);
    point3d <double> C(0, 5, 5);
    point3d <double> D(5, 0, 0);
    point3d <double> E(1, 1, 1);

    std::vector <double> len{1, 2, 3};
    point3d <double> F(1, 5, 3);
    point3d <double> G(6, 2, 4);
    
    const double epsilon = 1e-6;
    EXPECT_DOUBLE_EQ(weighed_euclidean(A, B, link_length), sqrt(125));
    EXPECT_DOUBLE_EQ(weighed_euclidean(A, C, link_length), 5 * sqrt(15));
    EXPECT_DOUBLE_EQ(weighed_euclidean(C, D, link_length), 5 * sqrt(30));
    EXPECT_DOUBLE_EQ(weighed_euclidean(A, E, link_length), sqrt(30));
    EXPECT_DOUBLE_EQ(weighed_euclidean(D, E, link_length), sqrt(255));
    EXPECT_DOUBLE_EQ(weighed_euclidean(F, G, len), 3 * sqrt(22));
    EXPECT_DOUBLE_EQ(weighed_euclidean(F, F, len), 0);
    EXPECT_DOUBLE_EQ(weighed_euclidean(G, G, len), 0);
}

// exceptions
TEST (weighedEuclidean, exceptions) {
    std::vector <double> link_length{5, 5};
    point3d <double> A(0, 0, 0);
    point2d <double> B(0, 5);
    point3d <double> C(0, 5, 5);
    point3d <double> D(5, 0, 0);

    EXPECT_THROW(weighed_euclidean(A, B, link_length), std::logic_error);
    EXPECT_THROW(weighed_euclidean(C, D, link_length), std::logic_error);
}

/*************************
on_line
*************************/

// test the standard 2D case (we´re not even using function in 3D)
TEST (online, 2D) {
    const double eps = 1e-10;
    std::pair <point2d<double>, point2d<double>> standard{point2d <double>{0, 0}, point2d <double>{5, 5}};
    point2d <double> A(2, 2);
    point2d <double> B(2-0.7*eps, 2+0.7*eps);
    point2d <double> C(1.999, 2.001);
    point2d <double> D(3, 2);
    point2d <double> E(-1,-1);
    point2d <double> F(6, 6);
    point2d <double> G(0, 0);

    EXPECT_TRUE(on_line<double>(A, standard));
    EXPECT_TRUE(on_line<double>(B, standard));
    EXPECT_FALSE(on_line<double>(C, standard));
    EXPECT_FALSE(on_line<double>(D, standard));
    EXPECT_FALSE(on_line<double>(E, standard));
    EXPECT_FALSE(on_line<double>(F, standard));
    EXPECT_TRUE(on_line<double>(G, standard));

    std::pair <point2d <double>, point2d <double>> upward{point2d <double>{0, 0}, point2d <double>{0, -1}};
    point2d <double> H(0, -0.2);
    point2d <double> I(0, 0.2);
    point2d <double> J(0, 0);
    point2d <double> K(0, -1);
    point2d <double> L(0.5*eps, -0.5*eps);
    point2d <double> M(0, -1+eps);
    point2d <double> N(-2*eps, 1);

    EXPECT_TRUE(on_line<double>(H, upward));
    EXPECT_FALSE(on_line<double>(I, upward));
    EXPECT_TRUE(on_line<double>(J, upward));
    EXPECT_TRUE(on_line<double>(K, upward));
    EXPECT_TRUE(on_line<double>(L, upward));
    EXPECT_TRUE(on_line<double>(M, upward));
    EXPECT_FALSE(on_line<double>(N, upward));
}

/*************************
lines_intersect2d
*************************/

// normal cases (nonparallel)
TEST (intersect2D, normal) {
    const double eps = 1e-6;

    // {(1, 6), (2, 5)}, {(1/3, 8), (2/3, 7)}
    point2d <double> A1(1, 6), A2(2, 5);
    std::pair <point2d <double>, point2d <double>> A{A1, A2};
    point2d <double> B1(double(1)/3, 8), B2(double(2)/3, 7);
    std::pair <point2d <double>, point2d <double>> B{B1, B2};
    std::pair <point2d <double>, point2d <double>> Bp{B2, A1};
    point2d <double> A3(1.0001, 5.9999);
    std::pair <point2d <double>, point2d <double>> Ap{A3, A2};

    EXPECT_FALSE(lines_intersect2d<double>(A, B));
    EXPECT_TRUE(lines_intersect2d<double>(A, Bp));
    EXPECT_FALSE(lines_intersect2d<double>(Ap, Bp));

    // {(1, 4), (4, 7)}, {(4, 4), (0, 12)}
    point2d <double> C1(1, 4), C2(4, 7);
    std::pair <point2d <double>, point2d <double>> C{C1, C2};
    point2d <double> D1(4, 4), D2(0, 12);
    std::pair <point2d <double>, point2d <double>> D{D1, D2};  
    point2d <double> C3(3.1, 6.1), C4(3, 6);
    std::pair <point2d <double>, point2d <double>> Cp{C3, C2};
    std::pair <point2d <double>, point2d <double>> Cq{C4, C2};

    EXPECT_TRUE(lines_intersect2d<double>(C, D));
    EXPECT_FALSE(lines_intersect2d<double>(Cp, D));
    EXPECT_TRUE(lines_intersect2d<double>(Cq, D));

    point2d <double> E1(1, 1), E2(4, 4);
    point2d <double> F1(1, 8), F2(3, 0);
    std::pair <point2d <double>, point2d <double>> E{E1, E2};
    std::pair <point2d <double>, point2d <double>> F{F1, F2};
    point2d <double> E3(2.41, 2.41);
    point2d <double> E4(2.4, 2.4);
    std::pair <point2d <double>, point2d <double>> Ep{E3, E2};
    std::pair <point2d <double>, point2d <double>> Eq{E4, E2};

    EXPECT_TRUE(lines_intersect2d<double>(E, F));
    EXPECT_FALSE(lines_intersect2d<double>(Ep, F));
    EXPECT_TRUE(lines_intersect2d<double>(Eq, F));

    point2d <double> G1(0, 0), G2(2, 2);
    point2d <double> H1(0, 2), H2(2, 0);
    std::pair <point2d <double>, point2d <double>> G{G1, G2};
    std::pair <point2d <double>, point2d <double>> H{H1, H2};
    point2d <double> G3(1.01, 1.01);
    point2d <double> G4(1, 1);
    std::pair <point2d <double>, point2d <double>> Gp{G3, G2};
    std::pair <point2d <double>, point2d <double>> Gq{G4, G2};

    EXPECT_TRUE(lines_intersect2d<double>(G, H));
    EXPECT_FALSE(lines_intersect2d<double>(Gp, H));
    EXPECT_TRUE(lines_intersect2d<double>(Gq, H));
}

// parallel lines should always give false
TEST (intersect2D, parallel) {
    // parallel
    point2d <double> A1(1, 6), A2(2, 5);
    point2d <double> B1(1, 5), B2(2, 4);
    std::pair <point2d <double>, point2d <double>> A{A1, A2};
    std::pair <point2d <double>, point2d <double>> B{B1, B2};
    EXPECT_FALSE(lines_intersect2d(A, B));

    const double eps = 1e-6;
    point2d <double> C1(1, 0), C2(2, 0);
    point2d <double> D1(1, 0.0001), D2(2, 0.0001), D3(1, 1e-12), D4(2, 1e-12);
    std::pair <point2d <double>, point2d <double>> C{C1, C2};
    std::pair <point2d <double>, point2d <double>> D{D1, D2};
    std::pair <point2d <double>, point2d <double>> Dp{D3, D4};

    EXPECT_FALSE(lines_intersect2d(C, D));
    EXPECT_TRUE(lines_intersect2d(C, Dp)); //overlapping
}

// if lines overlap give true
TEST (intersect2D, overlapping) {
    point2d <double> A1(0, 0), A2(0, 2);
    point2d <double> B1(0, -1), B2(0, 1), B3(0, 0), B4(0, -1e-6);
    std::pair <point2d <double>, point2d <double>> A{A1, A2};
    std::pair <point2d <double>, point2d <double>> B{B1, B2};
    std::pair <point2d <double>, point2d <double>> Bp{B1, B3};
    std::pair <point2d <double>, point2d <double>> Bq{B1, B4};
    
    EXPECT_TRUE(lines_intersect2d(A, B));
    EXPECT_TRUE(lines_intersect2d(A, Bp));
    EXPECT_FALSE(lines_intersect2d(A, Bq));

    point2d <double> C1(0, 0), C2(2, 4);
    point2d <double> D1(-1, -2), D2(1, 2);
    point2d <double> E1(-4, -2), E2(-0.5, -1);

    std::pair <point2d <double>, point2d <double>> C{C1, C2};
    std::pair <point2d <double>, point2d <double>> D{D1, D2};
    std::pair <point2d <double>, point2d <double>> E{E1, E2};
    EXPECT_TRUE(lines_intersect2d(C, D));
    EXPECT_FALSE(lines_intersect2d(C, E));
} 

/*************************
line_distance
*************************/

// line distance compared to online calculator https://keisan.casio.com/exec/system/1223531414
TEST (linedistance, nointersect) {
    const double eps = 1e-6;

    // line A (0, 1, 2) to (2, 3, 4)
    // line B (-4, 2, 3) to (3, -5, 1)
    point3d <double> A1(0, 1, 2), A2(2, 3, 4);
    point3d <double> B1(-4, 2, 3), B2(3, -5, 1);
    const double dist = 1.4385883441211;

    point3d <double> p0 = A1;
    point3d <double> r = A2 - A1;
    point3d <double> q0 = B1;
    point3d <double> s = B2 - B1;
    point3d <double> u = p0 - q0;

    EXPECT_NEAR(line_distance(r,s,u), dist, eps);
}

// distance for parallel lines
TEST (linedistance, parallel) {
    const double eps = 1e-6;

    // line C (0, 0, 0) to (2, 2, 2)
    // line D (2, 2, 1) to (0, 0, -1)
    point3d <double> C1(0, 0, 0), C2(2, 2, 2);
    point3d <double> D1(2, 2, 1), D2(0, 0, -1);
    const double dist = 0.81649658092773;

    point3d <double> p0 = C1;
    point3d <double> r = C2 - C1;
    point3d <double> q0 = D1;
    point3d <double> s = D2 - D1;
    point3d <double> u = p0 - q0;

    EXPECT_NEAR(line_distance(r,s,u), dist, eps);
}

// test intersecting lines example from http://mathforum.org/library/drmath/view/63719.html
TEST (linedistance, intersect) {
    const double eps = 1e-6;

    // line A (1, 0, 0) to (3, 3, 1)
    // line B (0, 5, 5) to (5, 6, 2)
    point3d <double> A1(1, 0, 0), A2(3, 3, 1);
    point3d <double> B1(0, 5, 5), B2(5, 6, 2);
    const double dist = 0;

    point3d <double> p0 = A1;
    point3d <double> r = A2 - A1;
    point3d <double> q0 = B1;
    point3d <double> s = B2 - B1;
    point3d <double> u = p0 - q0;

    EXPECT_NEAR(line_distance(r,s,u), dist, eps);
}

/*************************
signed_scaler
*************************/

TEST (signedscaler, all) {
    // line (1, 1, 1) to (3, 3, 3)
    point3d <double> A1(1, 1, 1), A2(3, 3, 3);  
    point3d <double> v = A2 - A1;
    point3d <double> u = A1 - A2;

    point3d <double> B(2, 2, 2);
    point3d <double> C(0, 0, 0);
    point3d <double> D(-1, -1, -1);
    point3d <double> E(3, 3, 3);

    EXPECT_DOUBLE_EQ(signed_scaler(B, v), 1);
    EXPECT_DOUBLE_EQ(signed_scaler(C, v), 0);
    EXPECT_DOUBLE_EQ(signed_scaler(D, v), -0.5);
    EXPECT_DOUBLE_EQ(signed_scaler(E, v), 1.5);

    EXPECT_DOUBLE_EQ(signed_scaler(B, u), -1);
    EXPECT_DOUBLE_EQ(signed_scaler(C, u), 0);
    EXPECT_DOUBLE_EQ(signed_scaler(D, u), 0.5);
    EXPECT_DOUBLE_EQ(signed_scaler(E, u), -1.5);
}

/*************************
line_intersects_box
*************************/

// standard intersections
TEST (intersectbox, standard) {
    std::vector <std::pair <double, double>> box{
        {0, 1}, {0, 1}, {0, 1}
    };

    // these intersect the box
    std::pair <point <double>, point <double>> linka {
        point<double>(std::vector <double>{-0.5, -0.5, -0.5}),
        point<double>(std::vector <double>{1.5, 1.5, 1.5})
    };
    std::pair <point <double>, point <double>> linkb {
        point<double>(std::vector <double>{0.5, 0.5, 1.2}),
        point<double>(std::vector <double>{1.2, 0.5, 0.5})
    };
    std::pair <point <double>, point <double>> linkc {
        point<double>(std::vector <double>{0.5, 0.5, -0.5}),
        point<double>(std::vector <double>{0.5, 0.5, 1.5})
    };
    std::pair <point <double>, point <double>> linkd {
        point<double>(std::vector <double>{1.5, 1.5, 1.5}),
        point<double>(std::vector <double>{-0.1, -0.2, -0.3})
    };

    // these dont intersect the box
    std::pair <point <double>, point <double>> linke {
        point<double>(std::vector <double>{0.5, 0.5, 1.5}),
        point<double>(std::vector <double>{1.5, 0.5, 1})
    };
    std::pair <point <double>, point <double>> linkf {
        point<double>(std::vector <double>{0.5, 0.5, 1.1}),
        point<double>(std::vector <double>{0.5, 0.5, 1.5})
    };
    std::pair <point <double>, point <double>> linkg {
        point<double>(std::vector <double>{1, 1, 1+1e-3}),
        point<double>(std::vector <double>{0, 1, 1+1e-3})
    };
    std::pair <point <double>, point <double>> linkh {
        point<double>(std::vector <double>{1, 1, -1.05}),
        point<double>(std::vector <double>{1, 1, -2})
    };

    EXPECT_TRUE(line_intersects_box(linka, box));
    EXPECT_TRUE(line_intersects_box(linkb, box));
    EXPECT_TRUE(line_intersects_box(linkc, box));
    EXPECT_TRUE(line_intersects_box(linkd, box));

    EXPECT_FALSE(line_intersects_box(linke, box));
    EXPECT_FALSE(line_intersects_box(linkf, box));
    EXPECT_FALSE(line_intersects_box(linkg, box));
    EXPECT_FALSE(line_intersects_box(linkh, box));
}

// test intersection where endpoints are inside box
TEST (intersectbox, inside) {
    std::vector <std::pair <double, double>> box{
        {0, 1}, {0, 1}, {0, 1}
    };

    // one endpoint inside box
    std::pair <point <double>, point <double>> linka {
        point<double>(std::vector <double>{1, 1, 1.5}),
        point<double>(std::vector <double>{0.5, 0.5, 0.5})
    };

    // both endpoints inside box
    std::pair <point <double>, point <double>> linkb {
        point<double>(std::vector <double>{0.5, 0.5, 0.5}),
        point<double>(std::vector <double>{0.2, 0.2, 0.2})
    };
    std::pair <point <double>, point <double>> linkc {
        point<double>(std::vector <double>{0.5, 0.5, 0.5}),
        point<double>(std::vector <double>{0, 0, 0})
    };

    EXPECT_TRUE(line_intersects_box(linka, box));
    EXPECT_TRUE(line_intersects_box(linkb, box));
    EXPECT_TRUE(line_intersects_box(linkc, box));
}

// test corner intersection
TEST (intersectbox, corners) {
    std::vector <std::pair <double, double>> box{
        {0, 1}, {0, 1}, {0, 1}
    };

    // corner of box
    std::pair <point <double>, point <double>> linka {
        point<double>(std::vector <double>{-0.5, -0.5, -0.5}),
        point<double>(std::vector <double>{1.5, 1.5, 1.5})
    };

    // edge of box
    std::pair <point <double>, point <double>> linkb {
        point<double>(std::vector <double>{-0.5, -0.5, 0.5}),
        point<double>(std::vector <double>{1.5, 1.5, 0.5})
    };
    std::pair <point <double>, point <double>> linkc {
        point<double>(std::vector <double>{-1, 1, 0.5}),
        point<double>(std::vector <double>{1, 1, 0.5})
    };

    EXPECT_TRUE(line_intersects_box(linka, box));
    EXPECT_TRUE(line_intersects_box(linkb, box));
    EXPECT_TRUE(line_intersects_box(linkc, box));
}

// test exceptions
TEST (intersectbox, exceptions) {
    std::vector <std::pair <double, double>> box{
        {0, 1}, {0, 1}
    };

    std::pair <point <double>, point <double>> linka {
        point<double>(std::vector <double>{-0.5, -0.5, 0.5}),
        point<double>(std::vector <double>{1.5, 1.5})
    };

    std::pair <point <double>, point <double>> linkb {
        point<double>(std::vector <double>{-0.5, -0.5, 0.5}),
        point<double>(std::vector <double>{1.5, 1.5, 1.5})
    };

    EXPECT_THROW(line_intersects_box(linka, box), std::logic_error);
    EXPECT_THROW(line_intersects_box(linkb, box), std::logic_error);
}

/*************************
lines_intersect_3d
*************************/

// retest 2D cases using 3D checking
TEST (intersect3D, planar) {
    {
        const double eps = 1e-6;

        // {(1, 6), (2, 5)}, {(1/3, 8), (2/3, 7)}
        point3d <double> A1(1, 6, 0), A2(2, 5, 0);
        std::pair <point3d <double>, point3d <double>> A{A1, A2};
        point3d <double> B1(double(1)/3, 8, 0), B2(double(2)/3, 7, 0);
        std::pair <point3d <double>, point3d <double>> B{B1, B2};
        std::pair <point3d <double>, point3d <double>> Bp{B2, A1};
        point3d <double> A3(1.0001, 5.9999, 0);
        std::pair <point3d <double>, point3d <double>> Ap{A3, A2};

        EXPECT_FALSE(lines_intersect3d<double>(A, B));
        EXPECT_TRUE(lines_intersect3d<double>(A, Bp));
        EXPECT_FALSE(lines_intersect3d<double>(Ap, Bp));

        // {(1, 4), (4, 7)}, {(4, 4), (0, 12)}
        point3d <double> C1(1, 4, 0), C2(4, 7, 0);
        std::pair <point3d <double>, point3d <double>> C{C1, C2};
        point3d <double> D1(4, 4, 0), D2(0, 12, 0);
        std::pair <point3d <double>, point3d <double>> D{D1, D2};  
        point3d <double> C3(3.1, 6.1, 0), C4(3, 6, 0);
        std::pair <point3d <double>, point3d <double>> Cp{C3, C2};
        std::pair <point3d <double>, point3d <double>> Cq{C4, C2};

        EXPECT_TRUE(lines_intersect3d<double>(C, D));
        EXPECT_FALSE(lines_intersect3d<double>(Cp, D));
        EXPECT_TRUE(lines_intersect3d<double>(Cq, D));

        point3d <double> E1(1, 1, 0), E2(4, 4, 0);
        point3d <double> F1(1, 8, 0), F2(3, 0, 0);
        std::pair <point3d <double>, point3d <double>> E{E1, E2};
        std::pair <point3d <double>, point3d <double>> F{F1, F2};
        point3d <double> E3(2.41, 2.41, 0);
        point3d <double> E4(2.4, 2.4, 0);
        std::pair <point3d <double>, point3d <double>> Ep{E3, E2};
        std::pair <point3d <double>, point3d <double>> Eq{E4, E2};

        EXPECT_TRUE(lines_intersect3d<double>(E, F));
        EXPECT_FALSE(lines_intersect3d<double>(Ep, F));
        EXPECT_TRUE(lines_intersect3d<double>(Eq, F));

        point3d <double> G1(0, 0, 0), G2(2, 2, 0);
        point3d <double> H1(0, 2, 0), H2(2, 0, 0);
        std::pair <point3d <double>, point3d <double>> G{G1, G2};
        std::pair <point3d <double>, point3d <double>> H{H1, H2};
        point3d <double> G3(1.01, 1.01, 0);
        point3d <double> G4(1, 1, 0);
        std::pair <point3d <double>, point3d <double>> Gp{G3, G2};
        std::pair <point3d <double>, point3d <double>> Gq{G4, G2};

        EXPECT_TRUE(lines_intersect3d<double>(G, H));
        EXPECT_FALSE(lines_intersect3d<double>(Gp, H));
        EXPECT_TRUE(lines_intersect3d<double>(Gq, H));
    }

    {
        // parallel
        point3d <double> A1(1, 6, 1), A2(2, 5, 1);
        point3d <double> B1(1, 5, 1), B2(2, 4, 1);
        std::pair <point3d <double>, point3d <double>> A{A1, A2};
        std::pair <point3d <double>, point3d <double>> B{B1, B2};
        EXPECT_FALSE(lines_intersect3d(A, B));

        const double eps = 1e-6;
        point3d <double> C1(1, 0, 1), C2(2, 0, 1);
        point3d <double> D1(1, 0.0001, 1), D2(2, 0.0001, 1), D3(1, 1e-12, 1), D4(2, 1e-12, 1);
        std::pair <point3d <double>, point3d <double>> C{C1, C2};
        std::pair <point3d <double>, point3d <double>> D{D1, D2};
        std::pair <point3d <double>, point3d <double>> Dp{D3, D4};

        EXPECT_FALSE(lines_intersect3d(C, D));
        EXPECT_TRUE(lines_intersect3d(C, Dp));
    }

    {    
        // overlapping
        point3d <double> A1(0, 0, 0), A2(0, 2, 0);
        point3d <double> B1(0, -1, 0), B2(0, 1, 0), B3(0, 0, 0), B4(0, -1e-6, 0);
        std::pair <point3d <double>, point3d <double>> A{A1, A2};
        std::pair <point3d <double>, point3d <double>> B{B1, B2};
        std::pair <point3d <double>, point3d <double>> Bp{B1, B3};
        std::pair <point3d <double>, point3d <double>> Bq{B1, B4};
        
        EXPECT_TRUE(lines_intersect3d(A, B));
        EXPECT_TRUE(lines_intersect3d(A, Bp));
        EXPECT_FALSE(lines_intersect3d(A, Bq));

        point3d <double> C1(0, 0, 0), C2(2, 4, 0);
        point3d <double> D1(-1, -2, 0), D2(1, 2, 0);
        point3d <double> E1(-4, -2, 0), E2(-0.5, -1, 0);

        std::pair <point3d <double>, point3d <double>> C{C1, C2};
        std::pair <point3d <double>, point3d <double>> D{D1, D2};
        std::pair <point3d <double>, point3d <double>> E{E1, E2};
        EXPECT_TRUE(lines_intersect3d(C, D));
        EXPECT_FALSE(lines_intersect3d(C, E));
    }
}

// tests from 3D line distance
TEST (intersect3D, old) {
    // line A (0, 1, 2) to (2, 3, 4)
    // line B (-4, 2, 3) to (3, -5, 1)
    point3d <double> A1(0, 1, 2), A2(2, 3, 4);
    point3d <double> B1(-4, 2, 3), B2(3, -5, 1);
    std::pair <point3d <double>, point3d <double>> A{A1, A2};
    std::pair <point3d <double>, point3d <double>> B{B1, B2};
    EXPECT_FALSE(lines_intersect3d(A, B));

    // line C (0, 0, 0) to (2, 2, 2)
    // line D (2, 2, 1) to (0, 0, -1)
    point3d <double> C1(0, 0, 0), C2(2, 2, 2);
    point3d <double> D1(2, 2, 1), D2(0, 0, -1);
    std::pair <point3d <double>, point3d <double>> C{C1, C2};
    std::pair <point3d <double>, point3d <double>> D{D1, D2};
    EXPECT_FALSE(lines_intersect3d(C, D));

    // line E (0, -1.5, -0.5) to (3, 3, 1)
    // line F (0, 5, 5) to (5, 6, 2)
    point3d <double> E1(0, -1.5, -0.5), E2(6, 7.5, 2.5);
    point3d <double> F1(0, 5, 5), F2(5, 6, 2);
    std::pair <point3d <double>, point3d <double>> E{E1, E2};
    std::pair <point3d <double>, point3d <double>> F{F1, F2};
    EXPECT_TRUE(lines_intersect3d(E, F));
}

// test some parallel 3D lines
TEST (intersect3D, parallel) {
    // line A (0, 0, 0) to (1, 1, 1)
    // line B (1, 0, 0) to (2, 1, 1)
    point3d <double> A1(0, 0, 0), A2(1, 1, 1);
    point3d <double> B1(1, 0, 0), B2(2, 1, 1);
    std::pair <point3d <double>, point3d <double>> A{A1, A2};
    std::pair <point3d <double>, point3d <double>> B{B1, B2};
    EXPECT_FALSE(lines_intersect3d(A, B));

    // line C (0, 0, 0) to (2, 2, 2)
    // line D (-1, -1, -1) to (1, 1, 1)
    point3d <double> C1(0, 0, 0), C2(2, 2, 2);
    point3d <double> D1(-1, -1, -1), D2(1, 1, 1);
    std::pair <point3d <double>, point3d <double>> C{C1, C2};
    std::pair <point3d <double>, point3d <double>> D{D1, D2};
    EXPECT_TRUE(lines_intersect3d(C, D));

    // line E (0, 0, 0) to (2, 2, 2)
    // line F (-1, -1, -1) to (0, 0, 0)
    point3d <double> E1(0, 0, 0), E2(2, 2, 2);
    point3d <double> F1(-1, -1, -1), F2(0, 0, 0), F3(-1e-4, -1e-4, -1e-4);
    std::pair <point3d <double>, point3d <double>> E{E1, E2};
    std::pair <point3d <double>, point3d <double>> F{F1, F2};
    std::pair <point3d <double>, point3d <double>> Fp{F1, F3};
    EXPECT_TRUE(lines_intersect3d(E, F));
    EXPECT_FALSE(lines_intersect3d(E, Fp));

    // line G (0, 0, 0) to (2, 2, 2)
    // line H (0.5, 0.5, 0.5) to (1, 1, 1)
    point3d <double> G1(0, 0, 0), G2(2, 2, 2);
    point3d <double> H1(0.5, 0.5, 0.5), H2(1, 1, 1);
    std::pair <point3d <double>, point3d <double>> G{G1, G2};
    std::pair <point3d <double>, point3d <double>> H{H1, H2};
    EXPECT_TRUE(lines_intersect3d(G, H));
}

/*************************
workspace2d::collides
*************************/

TEST (collides2D, threesegment) {
    point2d <double> base(0, 0);
    std::vector <double> link_lengths{2, 2, 4};
    std::vector <std::pair <double,double>> lims{
        {-M_PI, M_PI}, {-M_PI, M_PI}, {-M_PI, M_PI}
    };

    arm <double>* planar = new arm2d <double>(base, link_lengths, lims);
    std::vector <polygon <double>> obstacles;
    workspace2d <double> ws(obstacles, planar);

    point <double> A(
        std::vector <double>{0, 0, 0}
    );
    point <double> B(
        std::vector <double>{0, 0, M_PI - 1e-3}
    );
    point <double> C(
        std::vector <double>{0, 0, M_PI}
    );
    point <double> D(
        std::vector <double>{0, M_PI_2, 3*M_PI_4 - 1e-3}
    );
    point <double> E(
        std::vector <double>{0, M_PI_2, 3*M_PI_4 + 1e-3}
    );

    EXPECT_FALSE(ws.collides(A));
    EXPECT_FALSE(ws.collides(B));
    EXPECT_TRUE(ws.collides(C));
    EXPECT_FALSE(ws.collides(D));
    EXPECT_TRUE(ws.collides(E));
}

TEST (collides2D, foursegment) {
    point2d <double> base(0, 0);
    std::vector <double> link_lengths{5, 5, 5, 5};
    std::vector <std::pair <double,double>> lims{
        {-M_PI, M_PI}, {-M_PI, M_PI}, {-M_PI, M_PI}, {-M_PI, M_PI} 
    };

    arm <double>* planar = new arm2d <double>(base, link_lengths, lims);
    std::vector <polygon <double>> obstacles;
    workspace2d <double> ws(obstacles, planar);

    point <double> A(
        std::vector <double>{0, 0, 0, 0}
    );
    point <double> B(
        std::vector <double>{0, M_PI_2, M_PI_2, 0}
    );
    point <double> C(
        std::vector <double>{0, -M_PI_2, -M_PI_2-0.1, 0}
    );
    point <double> D(
        std::vector <double>{0, 2*M_PI/3, 2*M_PI/3, 0}
    );
    point <double> E(
        std::vector <double>{0, M_PI_2, -M_PI, 0}
    );
    point <double> F(
        std::vector <double>{0, -M_PI_2, -M_PI_2-1, 0}
    );

    EXPECT_FALSE(ws.collides(A));
    EXPECT_FALSE(ws.collides(B));
    EXPECT_FALSE(ws.collides(C));
    EXPECT_TRUE(ws.collides(D));
    EXPECT_TRUE(ws.collides(E));
    EXPECT_TRUE(ws.collides(F));

    // edge case : equilateral triangle
    point <double> G(
        std::vector <double>{2*M_PI/3, 2*M_PI/3, 2*M_PI/3 - 1e-3, 0}
    );
    point <double> H(
        std::vector <double>{2*M_PI/3, 2*M_PI/3, 2*M_PI/3 + 1e-3, 0}
    );

    EXPECT_FALSE(ws.collides(G));
    EXPECT_TRUE(ws.collides(H));
}

/*************************
workspace3d::collides
*************************/

TEST (collides3D, antropomorphic) {
    std::vector <box <double>> obstacles;

    point3d <double> base(0, 0, 0);
    std::vector <double> link_lengths{5, 5, 8};
    std::vector <std::pair <double,double>> lims{
        {-M_PI, M_PI}, {-M_PI, M_PI}, {-M_PI, M_PI}
    };

    arm <double>* antro = new antropomorphic_arm <double>(base, link_lengths, lims);
    workspace3d <double> ws(obstacles, antro);

    point <double> A(
        std::vector <double>{0, 0, 0}
    );
    point <double> B(
        std::vector <double>{0, 0, -M_PI}
    );
    point <double> C(
        std::vector <double>{M_PI, 0, -M_PI}
    );
    point <double> D(
        std::vector <double>{0, 0, -3*M_PI/4 - 1e-6}
    );
    point <double> E(
        std::vector <double>{M_PI_2, 0, -3*M_PI/4 - 1e-6}
    );
    point <double> F(
        std::vector <double>{M_PI_2, 0, -7*M_PI_2 / 5}
    );

    EXPECT_FALSE(ws.collides(A));
    EXPECT_TRUE(ws.collides(B));
    EXPECT_TRUE(ws.collides(C));
    EXPECT_TRUE(ws.collides(D));
    EXPECT_TRUE(ws.collides(E));
    EXPECT_FALSE(ws.collides(F));
}


int main (int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}