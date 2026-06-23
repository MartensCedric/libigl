#include "default_types.h"
#include <igl/winding_number_one_shot_2d.h>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/eigen/dense.h>

namespace nb = nanobind;
using namespace nb::literals;

namespace pyigl
{
  // Binding for winding_number_one_shot_2d over multiple query points
  auto winding_number_one_shot_2d(
    const nb::DRef<const Eigen::MatrixXN> &C,
    const nb::DRef<const Eigen::MatrixXN> &Q)
  {
    Eigen::VectorXN W;
    igl::winding_number_one_shot_2d(C, Q, W);
    return W;
  }
  // Binding for winding_number_one_shot_2d at a single query point
  Numeric winding_number_one_shot_2d_single(
    const nb::DRef<const Eigen::MatrixXN> &C,
    const nb::DRef<const Eigen::VectorXN> &q)
  {
    Eigen::MatrixXN Q = q.transpose();
    Eigen::VectorXN W;
    igl::winding_number_one_shot_2d(C, Q, W);
    return W(0);
  }
}

// Bind the wrapper to the Python module
void bind_winding_number_one_shot_2d(nb::module_ &m)
{
  m.def(
    "winding_number_one_shot_2d",
    &pyigl::winding_number_one_shot_2d,
    "C"_a,
    "Q"_a,
    R"(Generalized winding number for a 2D cubic Bézier spline via the One-Shot
    method [Martens & Bessmeltsev, "One-Shot Generalized Winding Numbers", CGF 2026].

    Computes the generalized winding number at each query point with respect to a
    collection of N directed cubic Bézier curves. Closed CCW splines return +1
    inside and 0 outside; open curves return fractional values in (0, 1).

    @param[in] C  4*N by 2 matrix of cubic Bézier control points; rows 4*i..4*i+3
                  are P0..P3 of the i-th curve
    @param[in] Q  #Q by 2 matrix of 2D query points
    @return Vector of generalized winding numbers, one per query point)");
  m.def(
    "winding_number_one_shot_2d",
    &pyigl::winding_number_one_shot_2d_single,
    "C"_a,
    "q"_a,
    R"(Generalized winding number for a 2D cubic Bézier spline via the One-Shot
    method, evaluated at a single query point.

    @param[in] C  4*N by 2 matrix of cubic Bézier control points; rows 4*i..4*i+3
                  are P0..P3 of the i-th curve
    @param[in] q  2-vector query point
    @return generalized winding number at q)");
}
