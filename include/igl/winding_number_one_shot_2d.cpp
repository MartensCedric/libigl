// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2026 Cedric Martens <cedric.martens@umontreal.ca>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "winding_number_one_shot_2d.h"
#include "ray_bezier_intersect_2d.h"
#include "PI.h"
#include "parallel_for.h"

#include <Eigen/Core>
#include <cmath>
#include <vector>

template <typename DerivedC, typename DerivedQ, typename DerivedW>
IGL_INLINE void igl::winding_number_one_shot_2d(
  const Eigen::MatrixBase<DerivedC> & C,
  const Eigen::MatrixBase<DerivedQ> & Q,
  Eigen::PlainObjectBase<DerivedW> & W)
{
  using Scalar = typename DerivedQ::Scalar;
  using RowVec2 = Eigen::Matrix<Scalar, 1, 2>;

  assert(C.cols() == 2 && "C must be 4*N by 2");
  assert(C.rows() % 4 == 0 && "C must have a multiple of 4 rows");
  assert(Q.cols() == 2 && "Q must be #Q by 2");

  const int N = static_cast<int>(C.rows()) / 4;
  const int nq = static_cast<int>(Q.rows());

  const Scalar TWO_PI = Scalar(2) * Scalar(igl::PI);

  // Per-curve AABB for early ray rejection (convex-hull property of Bézier).
  // bbox[i] = {xmin, xmax, ymin, ymax}
  struct BBox { Scalar xmin, xmax, ymin, ymax; };
  std::vector<BBox> bboxes(N);
  for (int i = 0; i < N; ++i)
  {
    auto ctrl = C.block(4*i, 0, 4, 2);
    bboxes[i] = {
      static_cast<Scalar>(ctrl.col(0).minCoeff()),
      static_cast<Scalar>(ctrl.col(0).maxCoeff()),
      static_cast<Scalar>(ctrl.col(1).minCoeff()),
      static_cast<Scalar>(ctrl.col(1).maxCoeff())
    };
  }

  W.resize(nq, 1);

  igl::parallel_for(nq, [&](const int qi)
  {
    const RowVec2 q(static_cast<Scalar>(Q(qi,0)), static_cast<Scalar>(Q(qi,1)));
    Scalar w = Scalar(0);

    for (int i = 0; i < N; ++i)
    {
      const BBox & bb = bboxes[i];
      const RowVec2 a = C.row(4*i).template cast<Scalar>();
      const RowVec2 b = C.row(4*i+3).template cast<Scalar>();

      // Evaluate control points as a concrete Matrix<Scalar,4,2> so the call to
      // ray_bezier_intersect_2d resolves to a pre-instantiated specialization.
      const Eigen::Matrix<Scalar, 4, 2> ctrl =
        C.block(4*i, 0, 4, 2).template cast<Scalar>();

      // Signed ray-intersection count (AABB early reject).
      int chi = 0;
      if (bb.xmax >= q(0) && bb.ymin <= q(1) && bb.ymax >= q(1))
        chi = igl::ray_bezier_intersect_2d(ctrl, q(1), q(0));

      // Endpoint arc fraction (one-shot formula).
      const RowVec2 ra = a - q;
      const RowVec2 rb = b - q;

      if (ra.squaredNorm() < Scalar(1e-14) || rb.squaredNorm() < Scalar(1e-14))
      {
        // Query coincides with an endpoint — full-circle limit.
        w += static_cast<Scalar>(chi);
        continue;
      }

      auto norm_angle = [TWO_PI](Scalar ang) -> Scalar
      {
        ang = std::fmod(ang, TWO_PI);
        if (ang < Scalar(0)) { ang += TWO_PI; }
        return ang;
      };

      const Scalar theta_A = norm_angle(std::atan2(ra(1), ra(0)));
      const Scalar theta_B = norm_angle(std::atan2(rb(1), rb(0)));

      Scalar delta = theta_B - theta_A;
      if (delta <= Scalar(0)) { delta += TWO_PI; }

      if (delta < Scalar(1e-10) || TWO_PI - delta < Scalar(1e-10))
      {
        // Endpoints collinear with query — no arc splitting.
        w += static_cast<Scalar>(chi);
        continue;
      }

      // Inner arc spans CCW from A to B (angular span delta).
      // The +x ray lies in the inner arc iff theta_A > theta_B.
      const int chi_outer = (theta_A > theta_B) ? chi - 1 : chi;
      w += static_cast<Scalar>(chi_outer) + delta / TWO_PI;
    }

    W(qi) = w;
  }, 1000);
}

#ifdef IGL_STATIC_LIBRARY
template void igl::winding_number_one_shot_2d<
  Eigen::MatrixXd, Eigen::MatrixXd, Eigen::VectorXd>(
    const Eigen::MatrixBase<Eigen::MatrixXd>&,
    const Eigen::MatrixBase<Eigen::MatrixXd>&,
    Eigen::PlainObjectBase<Eigen::VectorXd>&);
template void igl::winding_number_one_shot_2d<
  Eigen::MatrixXf, Eigen::MatrixXf, Eigen::VectorXf>(
    const Eigen::MatrixBase<Eigen::MatrixXf>&,
    const Eigen::MatrixBase<Eigen::MatrixXf>&,
    Eigen::PlainObjectBase<Eigen::VectorXf>&);
#endif
