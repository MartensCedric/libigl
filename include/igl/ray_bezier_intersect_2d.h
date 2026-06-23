// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2026 Cedric Martens <cedric.martens@umontreal.ca>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_RAY_BEZIER_INTERSECT_2D_H
#define IGL_RAY_BEZIER_INTERSECT_2D_H

#include "igl_inline.h"
#include <Eigen/Core>

namespace igl
{
  /// Signed intersection count of a horizontal half-line with a cubic Bézier curve.
  ///
  /// Counts intersections of { (x, y0) : x >= x_min } with the cubic Bézier
  /// curve whose 4 control points are the rows of C (in order P0, P1, P2, P3).
  /// Each crossing contributes +1 if the curve moves upward (tangent y > 0)
  /// and -1 if it moves downward. Tangential touches (zero y-component of
  /// tangent) are dropped.
  ///
  /// Uses Bézier clipping (recursive subdivision) for robustness and efficiency.
  /// See "One-Shot Generalized Winding Numbers" [Martens & Bessmeltsev, CGF 2026].
  ///
  /// @param[in]  C      4 by 2 matrix of cubic Bézier control points (rows P0..P3)
  /// @param[in]  y0     y-coordinate of the horizontal half-line
  /// @param[in]  x_min  left endpoint of the half-line; only crossings at x ≥ x_min
  ///                    are counted (default: −∞, i.e. the full horizontal line)
  /// @param[in]  tol    subdivision convergence tolerance (default: 1e-9)
  /// @return Signed crossing count: +1 per upward crossing, −1 per downward crossing
  ///
  /// \see CubicBezierCurve2D
  /// \see WindingNumberOneShot2DScene
  /// \see winding_number_one_shot_2d
  template <typename DerivedC>
  IGL_INLINE int ray_bezier_intersect_2d(
    const Eigen::MatrixBase<DerivedC> & C,
    typename DerivedC::Scalar y0,
    typename DerivedC::Scalar x_min = typename DerivedC::Scalar(-1e300),
    typename DerivedC::Scalar tol   = typename DerivedC::Scalar(1e-9));
}

#ifndef IGL_STATIC_LIBRARY
#  include "ray_bezier_intersect_2d.cpp"
#endif

#endif
