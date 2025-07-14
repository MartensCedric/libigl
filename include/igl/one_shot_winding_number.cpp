#include "one_shot_winding_number.h"
#include "bezier_clip.h"
#include "ray_box_intersect.h"
#include "parallel_for.h"
#include "PI.h"
#include <vector>
#include <cassert>
#include <random>

template <
    typename DerivedA,
    typename DerivedX,
    typename DerivedW>
IGL_INLINE void igl::one_shot_winding_number(
    const Eigen::MatrixBase<DerivedA> &A,
    const Eigen::MatrixBase<DerivedX> &X,
    Eigen::PlainObjectBase<DerivedW> &W)
{
  assert(A.rows() == X.rows() && "Area and chi matrix must have matching number of rows.");
  assert(A.cols() == X.cols() && "Area and chi matrix must have matching number of columns.");
  W.resize(A.rows());

  for (int i = 0; i < A.rows(); i++)
  {
    W(i) = A.row(i).dot(X.row(i).template cast<typename DerivedA::Scalar>());
  }
}

template <typename DerivedQ, typename DerivedC, typename DerivedR, typename DerivedW>
IGL_INLINE void igl::one_shot_winding_number_cubic_bezier(
    Eigen::MatrixBase<DerivedQ> &Q,
    const Eigen::MatrixBase<DerivedC> &C,
    Eigen::PlainObjectBase<DerivedW> &W)
{
  int output_size = Q.rows();
  assert(output_size >= 2 && "A ray needs at least two points");

  W.resize(output_size);
  using Scalar = typename DerivedQ::Scalar;
  Eigen::Vector2<Scalar> start_point = C.row(0);
  Eigen::Vector2<Scalar> end_point = C.row(3);

  Eigen::Vector2<Scalar> dir = Q.row(1) - Q.row(0);
  dir.normalize();

  Eigen::Matrix2<typename DerivedQ::Scalar> bounds;

  const auto &p0 = C.row(0);
  const auto &p1 = C.row(1);
  const auto &p2 = C.row(2);
  const auto &p3 = C.row(3);

  bounds << p0.cwiseMin(p1).cwiseMin(p2).cwiseMin(p3),
      p0.cwiseMax(p1).cwiseMax(p2).cwiseMax(p3);

  auto [ts_sq, normals] = igl::bezier_clip(Q.row(0).transpose(), dir, C, 1e-8);

  std::vector<int> sign(ts_sq.size());
  for (int k = 0; k < normals.rows(); k++)
  {
    bool same_dir = normals.row(k).dot(dir) > 0.0;
    sign[k] = same_dir ? 1 : -1;
  }

  Eigen::Vector3<Scalar> n(0.0, 0.0, 1.0);
  Eigen::Vector3<Scalar> dir_3d(dir(0), dir(1), 0.0);
  Eigen::Vector2<Scalar> first_q = Q.row(0);

  for (int i = 0; i < Q.rows(); ++i)
  {
    Eigen::Vector2<Scalar> q = Q.row(i);
    bool is_inside_box = true;
    if(bounds(0, 0) > q(0) || bounds(0, 1) < q(0) || bounds(1,0) > q(1) || bounds(1,1) < q(1))
		    is_inside_box = false;

    if (!is_inside_box)
    {
      // compute linearization
      Eigen::Vector2<Scalar> V0 = start_point - q;
      Eigen::Vector2<Scalar> V1 = end_point - q;
      W(i) = 0.5 * M_1_PI * std::acos(V0.normalize().dot(V1.normalize()));
    }
    else
    {

      Eigen::Vector2<Scalar> dir_to_start = (start_point - q).normalized();
      Eigen::Vector2<Scalar> dir_to_end = (end_point - q).normalized();

      Scalar theta_start = std::acos(dir_to_start.dot(dir));
      Scalar theta_end = std::acos(dir_to_end.dot(dir));
      Scalar cos_theta = dir_to_start.dot(dir_to_end);
      Scalar theta_radians = std::acos(cos_theta);

      Eigen::Vector3<Scalar> dir_to_start_3d = Eigen::Vector3<Scalar>(dir_to_start(0), dir_to_start(1), 0.0);
      Eigen::Vector3<Scalar> dir_to_end_3d = Eigen::Vector3<Scalar>(dir_to_end(0), dir_to_end(1), 0.0);

      if (dir_to_start_3d.cross(dir_3d).dot(n) > 0.0)
        theta_start = 2.0 * M_PI - theta_start;

      if (dir_to_end_3d.cross(dir_3d).dot(n) > 0.0)
        theta_end = 2.0 * M_PI - theta_end;

      if (dir_to_start_3d.cross(dir_to_end_3d).dot(n) < 0.0)
      {
        theta_radians = 2.0 * M_PI - theta_radians;
      }

      Scalar current_t_sq = (q - first_q).squaredNorm();
      int chi = 0;
      for (int k = 0; k < ts_sq.rows(); k++)
      {
        if (current_t_sq <= ts_sq(k))
        {
          chi += sign[k];
        }
      }

      assert(theta_start >= 0.0);
      assert(theta_end >= 0.0);
      assert(theta_radians >= 0.0);

      int other_chi = chi;
      if (theta_start < theta_end) // inside
      {

        other_chi++;

        W(i) = (chi * (2.0 * M_PI - theta_radians) + other_chi * theta_radians) / (2.0 * M_PI);
      }
      else
      {
        other_chi--;
        W(i) = ((chi * theta_radians) + other_chi * (2.0 * M_PI - theta_radians)) / (2.0 * M_PI);
      }
    }
  }
}

#ifdef IGL_STATIC_LIBRARY
template void igl::one_shot_winding_number<Eigen::Matrix<double, -1, 2, 0, -1, 2>, Eigen::Matrix<int, -1, 2, 0, -1, 2>, Eigen::Matrix<double, -1, 1, 0, -1, 1>>(const Eigen::MatrixBase<Eigen::Matrix<double, -1, 2, 0, -1, 2>>  &, const Eigen::MatrixBase<Eigen::Matrix<int, -1, 2, 0, -1, 2>>  &, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1>> &);
#endif