// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2026 Cedric Martens <cedric.martens@umontreal.ca>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "ray_bezier_intersect_2d.h"
#include "cubic.h"

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

template <typename DerivedC>
IGL_INLINE int igl::ray_bezier_intersect_2d(
  const Eigen::MatrixBase<DerivedC> & C,
  typename DerivedC::Scalar y0,
  typename DerivedC::Scalar x_min,
  typename DerivedC::Scalar tol)
{
  using Scalar = typename DerivedC::Scalar;
  using Vec2   = Eigen::Matrix<Scalar, 1, 2>;

  // Evaluate point (via igl::cubic) and tangent at t, sharing power terms.
  // The tangent is the analytic derivative of the cubic Bézier, computed
  // inline so the point and tangent reuse the same (1-t)/t power terms.
  auto eval_and_tangent = [&C](Scalar t, Vec2 & P, Vec2 & T) {
    igl::cubic(C, t, P);
    const Scalar s  = Scalar(1) - t;
    const Scalar s2 = s * s, t2 = t * t;
    T = Scalar(3) * (s2          * (C.row(1) - C.row(0))
                   + Scalar(2)*s*t * (C.row(2) - C.row(1))
                   + t2            * (C.row(3) - C.row(2)));
  };

  // Split a 4-point control polygon at t via de Casteljau.
  using CP4 = std::array<Vec2, 4>;
  auto split = [](const CP4 & P, Scalar t, CP4 & left, CP4 & right) {
    Vec2 Q0 = (Scalar(1)-t)*P[0] + t*P[1];
    Vec2 Q1 = (Scalar(1)-t)*P[1] + t*P[2];
    Vec2 Q2 = (Scalar(1)-t)*P[2] + t*P[3];
    Vec2 R0 = (Scalar(1)-t)*Q0   + t*Q1;
    Vec2 R1 = (Scalar(1)-t)*Q1   + t*Q2;
    Vec2 S0 = (Scalar(1)-t)*R0   + t*R1;
    left  = {P[0], Q0, R0, S0};
    right = {S0,   R1, Q2, P[3]};
  };

  // 1D Bézier clipping: returns false when there is no root in [0,1].
  // Sets lo/hi to the tightest sub-interval that could contain a root.
  static const double u[4] = {0.0, 1.0/3.0, 2.0/3.0, 1.0};
  auto bezier_clip_1d = [](const Scalar d[4], Scalar & lo, Scalar & hi) -> bool {
    lo =  Scalar(1);
    hi = Scalar(-1);
    for (int i = 0; i < 4; ++i)
    {
      for (int j = i+1; j < 4; ++j)
      {
        Scalar d0 = d[i], d1 = d[j];
        Scalar u0 = Scalar(u[i]), u1 = Scalar(u[j]);
        if (d0 == d1) {
          if (std::abs(d0) < Scalar(1e-14)) {
            lo = std::min(lo, u0); hi = std::max(hi, u1);
          }
          continue;
        }
        if ((d0 > Scalar(0)) == (d1 > Scalar(0)) &&
            d0 != Scalar(0) && d1 != Scalar(0)) continue;
        Scalar t_cross = u0 + (u1-u0) * (-d0) / (d1-d0);
        lo = std::min(lo, t_cross);
        hi = std::max(hi, t_cross);
      }
    }
    if (hi < lo) return false;
    lo = std::max(lo, Scalar(0));
    hi = std::min(hi, Scalar(1));
    return lo <= hi;
  };

  struct Hit { Scalar x; int sign; };
  struct Task { CP4 seg; Scalar t_lo, t_hi; int depth; };

  std::vector<Hit>  results;
  std::vector<Task> stack;
  stack.reserve(64);

  CP4 c0 = {Vec2(C.row(0)), Vec2(C.row(1)), Vec2(C.row(2)), Vec2(C.row(3))};
  stack.push_back({c0, Scalar(0), Scalar(1), 0});

  const int   max_depth = 50;
  const Scalar tol_sqrt  = std::sqrt(tol);

  while (!stack.empty())
  {
    auto [seg, t_lo, t_hi, depth] = stack.back();
    stack.pop_back();

    Scalar d[4];
    for (int i = 0; i < 4; ++i) d[i] = seg[i](1) - y0;

    Scalar max_abs_d = std::max({std::abs(d[0]), std::abs(d[1]),
                                 std::abs(d[2]), std::abs(d[3])});
    Scalar t_span = t_hi - t_lo;

    if (max_abs_d < tol_sqrt || t_span < tol || depth >= max_depth)
    {
      Vec2 pos, tang;
      eval_and_tangent(Scalar(0.5) * (t_lo + t_hi), pos, tang);
      if (pos(0) >= x_min && std::abs(tang(1)) > Scalar(1e-14))
        results.push_back({pos(0), tang(1) > Scalar(0) ? +1 : -1});
      continue;
    }

    {
      bool all_pos = d[0]>Scalar(0)&&d[1]>Scalar(0)&&d[2]>Scalar(0)&&d[3]>Scalar(0);
      bool all_neg = d[0]<Scalar(0)&&d[1]<Scalar(0)&&d[2]<Scalar(0)&&d[3]<Scalar(0);
      if (all_pos || all_neg) continue;
    }

    Scalar u_lo, u_hi;
    if (!bezier_clip_1d(d, u_lo, u_hi)) continue;

    Scalar new_span = u_hi - u_lo;

    if (new_span < Scalar(1e-10))
    {
      // Skip root sitting at u = 1: next segment will pick it up at u = 0.
      if (u_hi >= Scalar(1) - Scalar(1e-10)) continue;
      Vec2 pos, tang;
      eval_and_tangent(t_lo + Scalar(0.5)*(u_lo + u_hi)*t_span, pos, tang);
      if (std::abs(pos(1) - y0) > tol_sqrt) continue;
      if (pos(0) >= x_min && std::abs(tang(1)) > Scalar(1e-14))
        results.push_back({pos(0), tang(1) > Scalar(0) ? +1 : -1});
      continue;
    }

    if (new_span > Scalar(0.8) || depth >= max_depth / 2)
    {
      // Slow progress: bisect instead of clip.
      CP4 left, right;
      split(seg, Scalar(0.5), left, right);
      Scalar t_mid_g = t_lo + Scalar(0.5) * t_span;
      stack.push_back({right, t_mid_g, t_hi,    depth+1});
      stack.push_back({left,  t_lo,    t_mid_g, depth+1});
    }
    else
    {
      // Clip to [u_lo, u_hi].
      CP4 clipped;
      if (u_lo > Scalar(1e-10))
      {
        CP4 lft, rgt;
        split(seg, u_lo, lft, rgt);
        Scalar t2 = (u_hi - u_lo) / (Scalar(1) - u_lo);
        CP4 cl, cr;
        split(rgt, t2, cl, cr);
        clipped = cl;
      }
      else
      {
        CP4 cl, cr;
        split(seg, u_hi, cl, cr);
        clipped = cl;
      }
      Scalar new_t_lo = t_lo + u_lo * t_span;
      Scalar new_t_hi = t_lo + u_hi * t_span;
      stack.push_back({clipped, new_t_lo, new_t_hi, depth+1});
    }
  }

  std::sort(results.begin(), results.end(),
    [](const Hit & a, const Hit & b){ return a.x < b.x; });

  // Merge near-duplicate hits and drop cancelled ones.
  std::vector<Hit> unique;
  unique.reserve(results.size());
  for (const auto & r : results)
  {
    if (!unique.empty() && std::abs(r.x - unique.back().x) < tol * Scalar(10))
      unique.back().sign += r.sign;
    else
      unique.push_back(r);
  }

  int total = 0;
  for (const auto & h : unique)
    if (h.sign != 0) total += h.sign;
  return total;
}

#ifdef IGL_STATIC_LIBRARY
template int igl::ray_bezier_intersect_2d<Eigen::Matrix<double,4,2>>(
  const Eigen::MatrixBase<Eigen::Matrix<double,4,2>>&, double, double, double);
template int igl::ray_bezier_intersect_2d<Eigen::Matrix<float,4,2>>(
  const Eigen::MatrixBase<Eigen::Matrix<float,4,2>>&, float, float, float);
template int igl::ray_bezier_intersect_2d<Eigen::MatrixXd>(
  const Eigen::MatrixBase<Eigen::MatrixXd>&, double, double, double);
#endif
