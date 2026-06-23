#include <test_common.h>
#include <igl/winding_number_one_shot_2d.h>
#include <igl/PI.h>

#include <Eigen/Core>
#include <cmath>
#include <random>

namespace
{
  // 4-arc cubic Bézier approximation of the unit circle (CCW).
  // Returns a 16×2 control-point matrix (4 rows per arc).
  Eigen::MatrixXd unit_circle_bezier()
  {
    const double k = 0.5522847498307936; // 4/3 * tan(pi/8)
    Eigen::MatrixXd C(16, 2);
    C <<  1, 0,   1, k,   k, 1,   0, 1,
          0, 1,  -k, 1,  -1, k,  -1, 0,
         -1, 0,  -1,-k,  -k,-1,   0,-1,
          0,-1,   k,-1,   1,-k,   1, 0;
    return C;
  }

  // Single open quarter arc from (1,0) to (0,1).
  Eigen::MatrixXd quarter_arc_bezier()
  {
    const double k = 0.5522847498307936;
    Eigen::MatrixXd C(4, 2);
    C << 1, 0,  1, k,  k, 1,  0, 1;
    return C;
  }
}

TEST_CASE("winding_number_one_shot_2d: closed circle interior", "[igl]")
{
  Eigen::MatrixXd C = unit_circle_bezier();
  Eigen::MatrixXd Q(1, 2); Q << 0, 0;
  Eigen::VectorXd W;
  igl::winding_number_one_shot_2d(C, Q, W);
  REQUIRE(W(0) == Approx(1.0).margin(1e-3));
}

TEST_CASE("winding_number_one_shot_2d: closed circle exterior", "[igl]")
{
  Eigen::MatrixXd C = unit_circle_bezier();
  Eigen::MatrixXd Q(2, 2);
  Q << 5, 5,
       1.5, 0;
  Eigen::VectorXd W;
  igl::winding_number_one_shot_2d(C, Q, W);
  REQUIRE(W(0) == Approx(0.0).margin(1e-3));
  REQUIRE(W(1) == Approx(0.0).margin(1e-3));
}

TEST_CASE("winding_number_one_shot_2d: open arc gives fractional value", "[igl]")
{
  Eigen::MatrixXd C = quarter_arc_bezier();
  Eigen::MatrixXd Q(1, 2); Q << 0.3, 0.3;
  Eigen::VectorXd W;
  igl::winding_number_one_shot_2d(C, Q, W);
  REQUIRE(std::abs(W(0)) > 1e-6);
  REQUIRE(std::abs(W(0)) < 1.0);
}

TEST_CASE("winding_number_one_shot_2d: batch matches single-point", "[igl]")
{
  Eigen::MatrixXd C = unit_circle_bezier();
  Eigen::MatrixXd Q(5, 2);
  Q << 0,0,  0.5,0.5,  -0.5,0.5,  3,3,  0.9,0;

  Eigen::VectorXd W_batch;
  igl::winding_number_one_shot_2d(C, Q, W_batch);

  for (int i = 0; i < Q.rows(); ++i)
  {
    Eigen::MatrixXd qi = Q.row(i);
    Eigen::VectorXd wi;
    igl::winding_number_one_shot_2d(C, qi, wi);
    REQUIRE(W_batch(i) == Approx(wi(0)).margin(1e-12));
  }
}

TEST_CASE("winding_number_one_shot_2d: float scalar compiles and is accurate", "[igl]")
{
  const float k = 0.5522847498307936f;
  Eigen::MatrixXf C(16, 2);
  C <<  1,0,  1,k,  k,1,  0,1,
        0,1, -k,1, -1,k, -1,0,
       -1,0, -1,-k, -k,-1, 0,-1,
        0,-1,  k,-1,  1,-k, 1,0;

  Eigen::MatrixXf Q(1, 2); Q << 0.f, 0.f;
  Eigen::VectorXf W;
  igl::winding_number_one_shot_2d(C, Q, W);
  REQUIRE(std::abs(W(0) - 1.f) < 1e-2f);
}

TEST_CASE("winding_number_one_shot_2d: query on endpoint is finite", "[igl]")
{
  Eigen::MatrixXd C = unit_circle_bezier();
  // Query exactly on the first control point (P0 of first arc = (1,0))
  Eigen::MatrixXd Q(1, 2); Q << 1.0, 0.0;
  Eigen::VectorXd W;
  igl::winding_number_one_shot_2d(C, Q, W);
  REQUIRE(std::isfinite(W(0)));
}

TEST_CASE("winding_number_one_shot_2d: 100 random arcs consistency", "[igl]")
{
  // Checks that circle winding numbers converge for many random query points.
  Eigen::MatrixXd C = unit_circle_bezier();

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> coord(-2.0, 2.0);

  Eigen::MatrixXd Q(100, 2);
  for (int i = 0; i < 100; ++i)
    Q.row(i) << coord(rng), coord(rng);

  Eigen::VectorXd W;
  igl::winding_number_one_shot_2d(C, Q, W);

  for (int i = 0; i < 100; ++i)
  {
    double r = Q.row(i).norm();
    if (r < 0.9)
      REQUIRE(W(i) == Approx(1.0).margin(1e-2));
    else if (r > 1.1)
      REQUIRE(W(i) == Approx(0.0).margin(1e-2));
  }
}

// ---------------------------------------------------------------------------
// Benchmark (hidden in normal runs via IGL_DEBUG_OFF == "[!hide]")
// ---------------------------------------------------------------------------
TEST_CASE("winding_number_one_shot_2d: benchmark", "[igl]" IGL_DEBUG_OFF)
{
  const double k = 0.5522847498307936;
  auto make_circle = [&](double cx, double cy, double r) -> Eigen::MatrixXd
  {
    Eigen::MatrixXd C(16, 2);
    C <<  cx+r,cy,    cx+r,cy+r*k, cx+r*k,cy+r, cx,cy+r,
          cx,cy+r,    cx-r*k,cy+r, cx-r,cy+r*k, cx-r,cy,
          cx-r,cy,    cx-r,cy-r*k, cx-r*k,cy-r, cx,cy-r,
          cx,cy-r,    cx+r*k,cy-r, cx+r,cy-r*k, cx+r,cy;
    return C;
  };

  std::mt19937 rng(1234);
  std::uniform_real_distribution<double> pos(-4.0, 4.0);
  std::uniform_real_distribution<double> rad(0.1, 0.5);

  // Build 250 circles (1000 arcs total)
  Eigen::MatrixXd C(250*16, 2);
  for (int i = 0; i < 250; ++i)
    C.block(i*16, 0, 16, 2) = make_circle(pos(rng), pos(rng), rad(rng));

  const int nq = 10000;
  std::uniform_real_distribution<double> qcoord(-5.0, 5.0);
  Eigen::MatrixXd Q(nq, 2);
  for (int i = 0; i < nq; ++i)
    Q.row(i) << qcoord(rng), qcoord(rng);

  BENCHMARK("one-shot 1000 arcs x 10000 queries [parallel batch]") {
    Eigen::VectorXd W;
    igl::winding_number_one_shot_2d(C, Q, W);
    return W.sum();
  };

  BENCHMARK("one-shot 1000 arcs x 10000 queries [serial loop]") {
    double sum = 0;
    for (int q = 0; q < nq; ++q) {
      Eigen::MatrixXd qi = Q.row(q);
      Eigen::VectorXd wi;
      igl::winding_number_one_shot_2d(C, qi, wi);
      sum += wi(0);
    }
    return sum;
  };
}
