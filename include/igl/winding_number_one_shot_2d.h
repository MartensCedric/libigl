// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2026 Cedric Martens <cedric.martens@umontreal.ca>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_WINDING_NUMBER_ONE_SHOT_2D_H
#define IGL_WINDING_NUMBER_ONE_SHOT_2D_H

#include "igl_inline.h"
#include <Eigen/Core>

namespace igl
{
  /// Generalized winding number for a 2D cubic Bézier spline via the One-Shot
  /// method [Martens & Bessmeltsev, "One-Shot Generalized Winding Numbers",
  /// CGF 2026, https://onlinelibrary.wiley.com/doi/full/10.1111/cgf.70194].
  ///
  /// Computes the generalized winding number at each query point with respect
  /// to a collection of N directed cubic Bézier curves.  A single horizontal
  /// ray is shot per curve per query; the signed crossing count is combined
  /// with the arc fraction subtended by the curve's endpoints to give the
  /// exact GWN.  Closed CCW splines return +1 inside and 0 outside; open
  /// curves return fractional values in (0, 1).
  ///
  /// Queries are parallelized via igl::parallel_for.  Per-curve axis-aligned
  /// bounding boxes (convex-hull property of Bézier) are used to skip the
  /// Bézier-clipping step for curves that cannot intersect the ray.
  ///
  /// \code{.cpp}
  /// // 4-Bézier approximation of a unit circle (CCW), query inside/outside
  /// const double k = 0.5522847498307936;
  /// Eigen::MatrixXd C(16, 2);
  /// C <<  1, 0,   1, k,   k, 1,   0, 1,   // arc Q1: (1,0)→(0,1)
  ///       0, 1,  -k, 1,  -1, k,  -1, 0,   // arc Q2: (0,1)→(-1,0)
  ///      -1, 0,  -1,-k,  -k,-1,   0,-1,   // arc Q3: (-1,0)→(0,-1)
  ///       0,-1,   k,-1,   1,-k,   1, 0;   // arc Q4: (0,-1)→(1,0)
  /// Eigen::MatrixXd Q(2, 2);
  /// Q << 0.0, 0.0,   // inside  → W ≈ 1
  ///      5.0, 5.0;   // outside → W ≈ 0
  /// Eigen::VectorXd W;
  /// igl::winding_number_one_shot_2d(C, Q, W);
  /// \endcode
  ///
  /// @param[in]  C  4*N by 2 matrix of cubic Bézier control points.
  ///               Rows 4*i … 4*i+3 are P0..P3 of the i-th curve.
  /// @param[in]  Q  #Q by 2 matrix of 2D query points.
  /// @param[out] W  #Q vector of generalized winding numbers.
  ///
  /// \see ray_bezier_intersect_2d
  /// \see igl::predicates::spline_winding_number
  template <typename DerivedC, typename DerivedQ, typename DerivedW>
  IGL_INLINE void winding_number_one_shot_2d(
    const Eigen::MatrixBase<DerivedC> & C,
    const Eigen::MatrixBase<DerivedQ> & Q,
    Eigen::PlainObjectBase<DerivedW> & W);
}

#ifndef IGL_STATIC_LIBRARY
#  include "winding_number_one_shot_2d.cpp"
#endif

#endif
