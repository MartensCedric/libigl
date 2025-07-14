#ifndef BEZIER_CLIP_H
#define BEZIER_CLIP_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <vector>
namespace igl
{

  template<typename DerivedQ, typename DerivedR, typename DerivedC, typename DerivedN, typename Scalar>
  std::pair<std::vector<double>, Eigen::MatrixBase<DerivedN>> bezier_clip(const DerivedQ& origin, const DerivedR& direction, const DerivedC& control_points, Scalar tolerance);


  template<typename DerivedQ, typename DerivedC, typename Scalar>
  std::pair<std::vector<typename DerivedQ::Scalar>, Eigen::MatrixBase<DerivedQ>> igl::bezier_clip(const DerivedQ& origin, const DerivedQ& direction, const DerivedC& control_points, Scalar t_min, Scalar t_max, Scalar tolerance);
}

#ifndef IGL_STATIC_LIBRARY
  #include "bezier_clip.cpp"
#endif

#endif