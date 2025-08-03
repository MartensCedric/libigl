#include "bezier_clip.h"
#include "bezier.h"
#include "ray_box_intersect.h"
#include "bezier_subdivide.h"
#include "PlainMatrix.h"
#include <cassert>

#include <iostream>

 template<typename DerivedQ, typename DerivedR, typename DerivedC, typename Scalar, typename DerivedT, typename DerivedN>
  IGL_INLINE void igl::bezier_clip(
    const Eigen::MatrixBase<DerivedQ>& origin, 
    const Eigen::MatrixBase<DerivedR>& direction, 
    const Eigen::MatrixBase<DerivedC>& control_points, 
    Scalar tolerance, 
    Eigen::PlainObjectBase<DerivedT>& t_sq, 
    Eigen::MatrixBase<DerivedN>& normals)
{

  //  static_assert(DerivedQ::RowsAtCompileTime == 1 && DerivedQ::ColsAtCompileTime == 2,
  //               "Origin must be a row vector of size 1x2");
  //   static_assert(DerivedR::RowsAtCompileTime == 1 && DerivedR::ColsAtCompileTime == 2,
  //               "Direction must be a row vector of size 1x2");

  igl::bezier_clip(origin, direction, control_points, 0.0, 1.0, tolerance, t_sq, normals);
}


  template<typename DerivedQ, typename DerivedR, typename DerivedC, typename Scalar, typename DerivedT, typename DerivedN>
  IGL_INLINE void igl::bezier_clip(
    const Eigen::MatrixBase<DerivedQ>& origin, 
    const Eigen::MatrixBase<DerivedR>& direction, 
    const Eigen::MatrixBase<DerivedC>& control_points, 
    Scalar t_min, 
    Scalar t_max, 
    Scalar tolerance, 
    Eigen::PlainObjectBase<DerivedT>& t_sq, 
    Eigen::MatrixBase<DerivedN>& normals)
{
  static_assert(DerivedQ::RowsAtCompileTime == 1 && DerivedQ::ColsAtCompileTime == 2, 
              "Origin must be a row vector of size 1x2");
  static_assert(DerivedR::RowsAtCompileTime == 1 && DerivedR::ColsAtCompileTime == 2, 
              "Direction must be a row vector of size 1x2");

  Eigen::Matrix2<typename DerivedQ::Scalar> bounds;

  const auto &p0 = control_points.row(0);
  const auto &p1 = control_points.row(1);
  const auto &p2 = control_points.row(2);
  const auto &p3 = control_points.row(3);

  bounds << p0.cwiseMin(p1).cwiseMin(p2).cwiseMin(p3),
      p0.cwiseMax(p1).cwiseMax(p2).cwiseMax(p3);

  Eigen::Vector2<typename DerivedQ::Scalar> origin_col(origin(0), origin(1));
  Eigen::Vector2<typename DerivedR::Scalar> direction_col(direction(0), direction(1));
  
  if (!ray_box_intersect(origin_col, direction_col, bounds))
  {
    return;
  }

  // If the curve is small enough, consider it an intersection
  // but also make sure we're far enough from the endpoints points.
  Scalar tol_sq = tolerance * tolerance;

  if ((p3 - p0).squaredNorm() <= tol_sq)
  {
    Scalar t = (t_min + t_max) / 2.0;

    Eigen::MatrixXd point;
    Eigen::MatrixXd tangent;

    igl::bezier(control_points, t, point, tangent);

    Eigen::Vector2<Scalar> normal(tangent(1), -tangent(0)); // Perpendicular to the tangent
    normals.resize(1, 2);
    normals.row(0) = normal;

    t_sq.resize(1);
    t_sq(0) = (point - origin).norm();

    return;
  }

  // Subdivide the curve
  Scalar midT = (t_min + t_max) / 2.0;
  DerivedC control_points1;
  DerivedC control_points2;
  igl::bezier_subdivide(control_points, midT, control_points1, control_points2);

  // Recursively clip the left and right sub-curves
  Eigen::VectorX<Scalar> t_left, t_right;
  Eigen::Matrix<Scalar, Eigen::Dynamic, 2> n_left, n_right;

  bezier_clip(origin, direction, control_points1, t_min, midT, tolerance, t_left, n_left);
  bezier_clip(origin, direction, control_points2, midT, t_max, tolerance, t_right, n_right);

  // Combine results
  t_sq.resize(t_left.size() + t_right.size());
  t_sq << t_left, t_right;

  normals.resize(n_left.rows() + n_right.rows(), 2);
  if (n_left.rows() > 0 && n_right.rows() > 0)
  {
    normals << n_left, n_right;
  }
  else if (n_left.rows() > 0)
  {
    normals = n_left;
  }
  else if (n_right.rows() > 0)
  {
    normals = n_right;
  }
  else
  {
    normals = Eigen::Matrix<Scalar, Eigen::Dynamic, 2>();
  }
}

#ifdef IGL_STATIC_LIBRARY

template void igl::bezier_clip<
    Eigen::Matrix<double, 1, 2>, 
    Eigen::Matrix<double, 1, 2>, 
    Eigen::Matrix<double, -1, -1>, 
    double, 
    Eigen::Matrix<double, -1, 1>, 
    Eigen::Matrix<double, -1, -1>
    >(
      const Eigen::MatrixBase<Eigen::Matrix<double, 1, 2>>&, 
      const Eigen::MatrixBase<Eigen::Matrix<double, 1, 2>>&, 
      const Eigen::MatrixBase<Eigen::Matrix<double, -1,-1> >&, 
      double, 
      Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1>>&, 
      Eigen::MatrixBase<Eigen::Matrix<double, -1, -1>>&);

template void igl::bezier_clip<
    Eigen::Matrix<double, 1, 2>, 
    Eigen::Matrix<double, 1, 2>, 
    Eigen::Matrix<double, -1, -1>, 
    double,
    Eigen::Matrix<double, -1, 1>, 
    Eigen::Matrix<double, -1, -1>
    >(
      const Eigen::MatrixBase<Eigen::Matrix<double, 1, 2>>&, 
      const Eigen::MatrixBase<Eigen::Matrix<double, 1, 2>>&, 
      const Eigen::MatrixBase<Eigen::Matrix<double, -1,-1> >&, 
      double,
      double,
      double,
      Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1>>&, 
      Eigen::MatrixBase<Eigen::Matrix<double, -1, -1>>&);
#endif
