#include <test_common.h>
#include <igl/one_shot_winding_number.h>
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
  Eigen::MatrixXd V;

}
