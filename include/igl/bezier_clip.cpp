#include "bezier_clip.h"
#include "bezier.h"
#include "ray_box_intersect.h"
#include "bezier_subdivide.h"
#include "PlainMatrix.h"
#include <cassert>

template <typename DerivedQ, typename DerivedR, typename DerivedC, typename DerivedN, typename Scalar>
IGL_INLINE void igl::bezier_clip(const Eigen::MatrixBase<DerivedQ> &origin, const Eigen::MatrixBase<DerivedR> &direction, const Eigen::MatrixBase<DerivedC> &control_points, Scalar tolerance, Eigen::VectorX<Scalar> &t_sq, Eigen::MatrixBase<DerivedN> &normals)
{
  igl::bezier_clip(origin, direction, control_points, 0.0, 1.0, tolerance, t_sq, normals);
}

template <typename DerivedQ, typename DerivedR, typename DerivedC, typename DerivedN, typename Scalar>
IGL_INLINE void igl::bezier_clip(const Eigen::MatrixBase<DerivedQ> &origin, const Eigen::MatrixBase<DerivedR> &direction, const Eigen::MatrixBase<DerivedC> &control_points, Scalar t_min, Scalar t_max, Scalar tolerance, Eigen::VectorX<Scalar> &t_sq, Eigen::MatrixBase<DerivedN> &normals)
{
  Eigen::Matrix2<typename DerivedQ::Scalar> bounds;

  const auto &p0 = control_points.row(0);
  const auto &p1 = control_points.row(1);
  const auto &p2 = control_points.row(2);
  const auto &p3 = control_points.row(3);

  bounds << p0.cwiseMin(p1).cwiseMin(p2).cwiseMin(p3),
      p0.cwiseMax(p1).cwiseMax(p2).cwiseMax(p3);

  // Check if ray intersects bounding box

  if (!ray_box_intersect(origin, direction, bounds))
  {
    return;
  }

  // If the curve is small enough, consider it an intersection
  // but also make sure we're far enough from the endpoints points.
  Scalar tol_sq = tolerance * tolerance;

  if ((p3 - p0).squaredNorm() <= tol_sq)
  {
    Scalar t = (t_min + t_max) / 2.0;

    Eigen::Vector2<Scalar> point;
    Eigen::Vector2<Scalar> tangent;

    igl::bezier(control_points, t, point, tangent);

    Eigen::Vector2<Scalar> normal(tangent(1), -tangent(0)); // Perpendicular to the tangent
    normals.resize(1, 2);
    normals.row(0) = normal;

    t_sq.resize(1);
    t_sq.row(0) = (point - origin).norm();

    return;
  }

  // Subdivide the curve
  Scalar midT = (t_min + t_max) / 2.0;
  DerivedC control_points1;
  DerivedC control_points2;
  igl::bezier_subdivide(control_points, midT, control_points1, control_points2);

  // Recursively clip the left and right sub-curves
  Eigen::VectorX<Scalar> t_left, t_right;
  Eigen::MatrixX<DerivedN> n_left, n_right;

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
    normals = Eigen::MatrixX<typename DerivedQ::Scalar>();
  }
}

#ifdef IGL_STATIC_LIBRARY
template void igl::bezier_clip<
    Eigen::Matrix<double, 1, 2, 1, 1, 2>,
    Eigen::Matrix<double, 1, 2, 1, 1, 2>,
    Eigen::Matrix<double, -1, -1>,
    Eigen::Matrix<double, -1, -1>,
    double>(
    const Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2>> &,
    const Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2>> &,
    const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1>> &,
    double,
    Eigen::VectorXd &,
    Eigen::MatrixBase<Eigen::Matrix<double, -1, -1>> &);

template void igl::bezier_clip<
    Eigen::Matrix<double, 1, 2, 1, 1, 2>, 
    Eigen::Matrix<double, 1, 2, 1, 1, 2>,
    Eigen::Matrix<double, -1, -1>,
    Eigen::Matrix<double, -1, -1>,
    double>(
    const Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2>> &,
    const Eigen::MatrixBase<Eigen::Matrix<double, 1, 2, 1, 1, 2>> &,
    const Eigen::MatrixBase<Eigen::Matrix<double, -1, -1>> &,
    double,
    double,
    double,
    Eigen::VectorXd &,
    Eigen::MatrixBase<Eigen::Matrix<double, -1, -1>> &);
#endif
