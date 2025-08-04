#include <test_common.h>
#include <igl/one_shot_winding_number.h>
#include <igl/winding_number.h>
#include <igl/bezier_clip.h>
#include <igl/bezier.h>
#include <iostream>

TEST_CASE("one_shot_winding_number: Weighted Sum", "[igl]")
{
  Eigen::MatrixX2d A(5,2);
  Eigen::VectorXd leftA(5);
  Eigen::VectorXd rightA(5);

  leftA << 0.3, 0.2, 0.1, 0.9, 0.5;
  rightA << 0.7, 0.8, 0.9, 0.1, 0.5;

  A.col(0) = leftA;
  A.col(1) = rightA;

  Eigen::MatrixX2i X(5,2);
  Eigen::VectorXi leftX(5);
  Eigen::VectorXi rightX(5);

  leftX << 0, 0, 1, 2, -1;
  rightX << 1, -1, 0, 1, 0;

  X.col(0) = leftX;
  X.col(1) = rightX;

  Eigen::VectorXd gt(5);
  gt << 0.7, -0.8, 0.1, 1.9, -0.5;

  Eigen::VectorXd W(5);
  igl::one_shot_winding_number(A, X, W);
  test_common::assert_near(W, gt, 1e-12);

}


TEST_CASE("one_shot_winding_number: cubic bezier", "[igl]")
{
  Eigen::MatrixXd C(4,2);
  C.row(0) = Eigen::RowVector2d(0.0, 0.0);
  C.row(1) = Eigen::RowVector2d(1.0, 0.0);
  C.row(2) = Eigen::RowVector2d(1.0, 1.0);
  C.row(3) = Eigen::RowVector2d(1.0, 0.0);
  

  const auto &p0 = C.row(0);
  const auto &p1 = C.row(1);
  const auto &p2 = C.row(2);
  const auto &p3 = C.row(3);

  Eigen::MatrixXd bounds(2,2);
  bounds << p0.cwiseMin(p1).cwiseMin(p2).cwiseMin(p3),
      p0.cwiseMax(p1).cwiseMax(p2).cwiseMax(p3);

  Eigen::MatrixXd Q(10, 2);

  for(int i = 0; i < 10; i++)
  {
    Q.row(i) = Eigen::RowVector2d(static_cast<double>(i) - 3.0, 0.5);
  }

  Eigen::Matrix<double, Eigen::Dynamic, 1> ts_sq;
  Eigen::MatrixXd normals(0, 2);
  Eigen::RowVector2d Q_row = Q.row(0).head<2>();
  Eigen::RowVector2d dir = (Q.row(1) - Q.row(0)).normalized(); 
  igl::bezier_clip(Q_row, dir, C, 1e-8, ts_sq, normals);

  std::vector<int> sign(ts_sq.rows());
  for (int k = 0; k < normals.rows(); k++)
  {
    bool same_dir = normals.row(k).dot(dir) > 0.0;
    sign[k] = same_dir ? 1 : -1;
  }

  Eigen::MatrixXi S = Eigen::Map<Eigen::MatrixXi>(sign.data(), sign.size(), 1);
  Eigen::MatrixXd T_sq = Eigen::Map<Eigen::MatrixXd>(ts_sq.data(), ts_sq.size(), 1);

  Eigen::MatrixXd W(10, 1);
  Eigen::MatrixXd endpoints(2,2);
  endpoints.row(0) = p0;
  endpoints.row(1) = p3;

  igl::one_shot_winding_number(Q, endpoints, T_sq, S, bounds, true, W);
  
  Eigen::MatrixXd W_no_linear(10, 1);
  igl::one_shot_winding_number(Q, endpoints, T_sq, S, bounds, false, W_no_linear);

  test_common::assert_near(W, W_no_linear);


  const int num_samples = 100;
  Eigen::MatrixXd sampled_points(num_samples, 2);

  for (int i = 0; i < num_samples; i++) {
    double t = static_cast<double>(i) / (num_samples - 1);
    Eigen::MatrixXd value;
    igl::bezier(C, t, value);
    sampled_points.row(i) = value;
  }

  Eigen::MatrixXi E(num_samples - 1, 2);
  for (int i = 0; i < num_samples - 1; i++) {
    E(i, 0) = i;
    E(i, 1) = i + 1;
  }

  Eigen::MatrixXd W_direct(10, 1);
  for (int i = 0; i < 10; i++) {
    Eigen::Matrix<double, 1, 2> q_i =  Q.row(i).head<2>();
    W_direct(i) = igl::winding_number(sampled_points, E, q_i);
  }

  test_common::assert_near(W, W_direct);


  std::cout << "W: " << W << std::endl;
  std::cout << "W_no_linear: " << W_no_linear << std::endl;
  std::cout << "W_direct: " << W_direct << std::endl;

  test_common::assert_eq(false, true);
}
