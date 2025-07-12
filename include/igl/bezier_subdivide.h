#ifndef BEZIER_SUBDIVIDE_H
#define BEZIER_SUBDIVIDE_H
#include "igl_inline.h"
#include <Eigen/Core>
#include <vector>
namespace igl
{

template <typename DerivedV>
void bezier_subdivide( const Eigen::MatrixBase<DerivedV> & C, const typename DerivedV::Scalar t, Eigen::MatrixBase<DerivedV> & C1, Eigen::MatrixBase<DerivedV> & C2);
}
#ifndef IGL_STATIC_LIBRARY
  #include "bezier_subdivide.cpp"
#endif

#endif

