#ifndef BEZIER_CLIP_H
#define BEZIER_CLIP_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <vector>
namespace igl
{

  template<typename DerivedQ, typename DerivedR, typename DerivedC, typename DerivedN>
  std::pair<std::vector<double>, Eigen::MatrixBase<DerivedN>> bezier_clip(const DerivedQ& origin, const DerivedR& direction, const DerivedC& control_points, double tolerance);


  template<typename DerivedQ, typename DerivedR, typename DerivedC, typename DerivedN>
  std::pair<std::vector<double>, Eigen::MatrixBase<DerivedN>> bezier_clip(const DerivedQ& origin, const DerivedR& direction, const DerivedC& control_points, double t_min, double t_max, double tolerance);
}

#ifndef IGL_STATIC_LIBRARY
  #include "bezier_clip.cpp"
#endif

#endif