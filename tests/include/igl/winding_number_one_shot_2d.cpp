#include <test_common.h>
#include <igl/winding_number_one_shot_2d.h>
#include <igl/winding_number.h>
#include <igl/bezier.h>

#include <Eigen/Core>
#include <cmath>
#include <random>
#include <vector>

namespace
{
  // 4-arc cubic Bézier approximation of a CCW circle centered at (cx,cy).
  // Returns a 16×2 control-point matrix (4 rows per arc).
  Eigen::MatrixXd circle_bezier(double cx = 0, double cy = 0, double r = 1)
  {
    const double k = 0.5522847498307936; // 4/3 * tan(pi/8)
    Eigen::MatrixXd C(16, 2);
    C <<  1, 0,   1, k,   k, 1,   0, 1,
          0, 1,  -k, 1,  -1, k,  -1, 0,
         -1, 0,  -1,-k,  -k,-1,   0,-1,
          0,-1,   k,-1,   1,-k,   1, 0;
    C *= r;
    C.col(0).array() += cx;
    C.col(1).array() += cy;
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

  // Subdivide each cubic Bézier in C (4N×2) into n line segments and return
  // the polyline as an edge mesh (V,E) for igl::winding_number. Sampling is
  // done with igl::bezier's spline overload.
  void bezier_polyline(
    const Eigen::MatrixXd & C,
    int n,
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & E)
  {
    const int N = static_cast<int>(C.rows()) / 4;
    std::vector<Eigen::MatrixXd> spline(N);
    for (int i = 0; i < N; ++i)
    {
      spline[i] = C.block(4*i, 0, 4, 2);
    }
    const Eigen::VectorXd T = Eigen::VectorXd::LinSpaced(n+1, 0.0, 1.0);
    igl::bezier(spline, T, V);
    E.resize(N*n, 2);
    for (int i = 0; i < N; ++i)
    {
      for (int j = 0; j < n; ++j)
      {
        E.row(i*n+j) << i*(n+1)+j, i*(n+1)+j+1;
      }
    }
  }

  // Distance from q to the nearest polyline vertex; used to reject query
  // points too close to the curve, where the polyline ground truth itself
  // has non-negligible discretization error.
  double distance_to_polyline(const Eigen::MatrixXd & V, const Eigen::RowVector2d & q)
  {
    return (V.rowwise() - q).rowwise().norm().minCoeff();
  }

  // Compare one-shot winding numbers of the Bézier collection C against
  // igl::winding_number on its dense polyline subdivision, at all rows of Q
  // farther than min_dist from the curve.
  void require_matches_polyline(
    const Eigen::MatrixXd & C,
    const Eigen::MatrixXd & Q,
    double min_dist,
    double tol)
  {
    const int n = 2000; // segments per Bézier
    Eigen::MatrixXd V;
    Eigen::MatrixXi E;
    bezier_polyline(C, n, V, E);

    Eigen::VectorXd W;
    igl::winding_number_one_shot_2d(C, Q, W);
    Eigen::VectorXd W_gt;
    igl::winding_number(V, E, Q, W_gt);

    // Compare only at queries far enough from the curve that the polyline
    // ground truth itself is accurate there.
    std::vector<int> keep;
    for (int i = 0; i < Q.rows(); ++i)
    {
      if (distance_to_polyline(V, Q.row(i)) >= min_dist) { keep.push_back(i); }
    }
    // Guard against the distance filter silently rejecting everything.
    REQUIRE(int(keep.size()) > Q.rows() / 2);
    Eigen::VectorXd W_keep(keep.size());
    Eigen::VectorXd W_gt_keep(keep.size());
    for (int i = 0; i < int(keep.size()); ++i)
    {
      W_keep(i) = W(keep[i]);
      W_gt_keep(i) = W_gt(keep[i]);
    }
    test_common::assert_near(W_keep, W_gt_keep, tol);
  }
}

TEST_CASE("winding_number_one_shot_2d: closed circle matches subdivided polyline", "[igl]")
{
  const Eigen::MatrixXd C = circle_bezier();

  // Regular grid spanning interior, exterior, and near-boundary points.
  const int gn = 15;
  Eigen::MatrixXd Q(gn*gn, 2);
  for (int i = 0; i < gn; ++i)
  {
    for (int j = 0; j < gn; ++j)
    {
      Q.row(i*gn+j) << -2.0 + 4.0*i/(gn-1), -2.0 + 4.0*j/(gn-1);
    }
  }
  require_matches_polyline(C, Q, 0.1, 1e-4);

  // Closed CCW curve: sanity-check the classic integer values directly.
  Eigen::MatrixXd Q2(2, 2);
  Q2 << 0, 0,
        5, 5;
  Eigen::VectorXd W;
  igl::winding_number_one_shot_2d(C, Q2, W);
  Eigen::VectorXd W_expected(2);
  W_expected << 1.0, 0.0;
  test_common::assert_near(W, W_expected, 1e-6);
}

TEST_CASE("winding_number_one_shot_2d: open arc matches subdivided polyline", "[igl]")
{
  const Eigen::MatrixXd C = quarter_arc_bezier();

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> coord(-1.5, 1.5);
  Eigen::MatrixXd Q(100, 2);
  for (int i = 0; i < Q.rows(); ++i)
  {
    Q.row(i) << coord(rng), coord(rng);
  }
  require_matches_polyline(C, Q, 0.1, 1e-4);

  // Open curves must give strictly fractional values.
  Eigen::MatrixXd Q2(1, 2);
  Q2 << 0.3, 0.3;
  Eigen::VectorXd W;
  igl::winding_number_one_shot_2d(C, Q2, W);
  REQUIRE(std::abs(W(0)) > 1e-6);
  REQUIRE(std::abs(W(0)) < 1.0);
}

TEST_CASE("winding_number_one_shot_2d: random open cubics match subdivided polyline", "[igl]")
{
  std::mt19937 rng(1234);
  std::uniform_real_distribution<double> coord(-1.0, 1.0);

  // A scene of 5 random (generally wild) open cubic Bézier curves.
  Eigen::MatrixXd C(5*4, 2);
  for (int i = 0; i < C.rows(); ++i)
  {
    C.row(i) << coord(rng), coord(rng);
  }

  std::uniform_real_distribution<double> qcoord(-2.0, 2.0);
  Eigen::MatrixXd Q(100, 2);
  for (int i = 0; i < Q.rows(); ++i)
  {
    Q.row(i) << qcoord(rng), qcoord(rng);
  }
  require_matches_polyline(C, Q, 0.15, 1e-3);
}

TEST_CASE("winding_number_one_shot_2d: overlapping circles match subdivided polyline", "[igl]")
{
  // Two overlapping CCW circles: winding number 2 in the lens, 1 in each
  // single-covered region, 0 outside.
  Eigen::MatrixXd C(32, 2);
  C.topRows(16) = circle_bezier(-0.5, 0.0, 1.0);
  C.bottomRows(16) = circle_bezier(0.5, 0.0, 1.0);

  Eigen::MatrixXd Q(4, 2);
  Q <<  0.0, 0.0,   // lens: w = 2
       -1.0, 0.0,   // left circle only: w = 1
        1.0, 0.0,   // right circle only: w = 1
        0.0, 3.0;   // outside: w = 0
  require_matches_polyline(C, Q, 0.1, 1e-4);

  Eigen::VectorXd W;
  igl::winding_number_one_shot_2d(C, Q, W);
  Eigen::VectorXd W_expected(4);
  W_expected << 2.0, 1.0, 1.0, 0.0;
  test_common::assert_near(W, W_expected, 1e-6);
}

TEST_CASE("winding_number_one_shot_2d: batch matches single-point", "[igl]")
{
  const Eigen::MatrixXd C = circle_bezier();
  Eigen::MatrixXd Q(5, 2);
  Q << 0,0,  0.5,0.5,  -0.5,0.5,  3,3,  0.9,0;

  Eigen::VectorXd W_batch;
  igl::winding_number_one_shot_2d(C, Q, W_batch);

  Eigen::VectorXd W_single(Q.rows());
  for (int i = 0; i < Q.rows(); ++i)
  {
    Eigen::MatrixXd qi = Q.row(i);
    Eigen::VectorXd wi;
    igl::winding_number_one_shot_2d(C, qi, wi);
    W_single(i) = wi(0);
  }
  test_common::assert_near(W_batch, W_single, 1e-12);
}

TEST_CASE("winding_number_one_shot_2d: float scalar compiles and is accurate", "[igl]")
{
  const Eigen::MatrixXf C = circle_bezier().cast<float>();
  Eigen::MatrixXf Q(1, 2);
  Q << 0.f, 0.f;
  Eigen::VectorXf W;
  igl::winding_number_one_shot_2d(C, Q, W);
  REQUIRE(std::abs(W(0) - 1.f) < 1e-2f);
}

TEST_CASE("winding_number_one_shot_2d: query on endpoint is finite", "[igl]")
{
  const Eigen::MatrixXd C = circle_bezier();
  // Query exactly on the first control point (P0 of first arc = (1,0))
  Eigen::MatrixXd Q(1, 2);
  Q << 1.0, 0.0;
  Eigen::VectorXd W;
  igl::winding_number_one_shot_2d(C, Q, W);
  REQUIRE(std::isfinite(W(0)));
}

// ---------------------------------------------------------------------------
// Benchmark (hidden in normal runs via IGL_DEBUG_OFF == "[!hide]")
// ---------------------------------------------------------------------------
TEST_CASE("winding_number_one_shot_2d: benchmark", "[igl]" IGL_DEBUG_OFF)
{
  std::mt19937 rng(1234);
  std::uniform_real_distribution<double> pos(-4.0, 4.0);
  std::uniform_real_distribution<double> rad(0.1, 0.5);

  // Build 250 circles (1000 arcs total)
  Eigen::MatrixXd C(250*16, 2);
  for (int i = 0; i < 250; ++i)
  {
    C.block(i*16, 0, 16, 2) = circle_bezier(pos(rng), pos(rng), rad(rng));
  }

  const int nq = 10000;
  std::uniform_real_distribution<double> qcoord(-5.0, 5.0);
  Eigen::MatrixXd Q(nq, 2);
  for (int i = 0; i < nq; ++i)
  {
    Q.row(i) << qcoord(rng), qcoord(rng);
  }

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
