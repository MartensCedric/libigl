#ifndef BEZIER_CLIP_H
#define BEZIER_CLIP_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <vector>
namespace igl
{
  template<typename DerivedQ, typename DerivedR, typename DerivedC, typename Scalar>
  void bezier_clip(
    const Eigen::MatrixBase<DerivedQ>& origin, 
    const Eigen::MatrixBase<DerivedR>& direction, 
    const Eigen::MatrixBase<DerivedC>& control_points, 
    Scalar tolerance, 
    Eigen::VectorX<Scalar>& t_sq, 
    Eigen::Matrix<Scalar, Eigen::Dynamic, 2>& normals);

  template<typename DerivedQ, typename DerivedR, typename DerivedC, typename Scalar>
  void bezier_clip(
    const Eigen::MatrixBase<DerivedQ>& origin, 
    const Eigen::MatrixBase<DerivedR>& direction, 
    const Eigen::MatrixBase<DerivedC>& control_points, 
    Scalar t_min, 
    Scalar t_max, 
    Scalar tolerance, 
    Eigen::VectorX<Scalar>& t_sq, 
    Eigen::Matrix<Scalar, Eigen::Dynamic, 2>& normals);
}

#ifndef IGL_STATIC_LIBRARY
  #include "bezier_clip.cpp"
#endif

#endif