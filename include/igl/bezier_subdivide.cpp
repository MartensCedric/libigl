#include "bezier_subdivide.h"
#include "igl_inline.h"
#include <Eigen/Core>
#include <vector>

template <typename DerivedV>
IGL_INLINE void bezier_subdivide( const Eigen::MatrixBase<DerivedV> & C, const typename DerivedV::Scalar t, Eigen::MatrixBase<DerivedV> & C1, Eigen::MatrixBase<DerivedV> & C2)
{
  typedef typename DerivedV::Scalar Scalar;
  const auto& p0 = C.row(0);
  const auto& p1 = C.row(1);
  const auto& p2 = C.row(2);
  const auto& p3 = C.row(3);

  Eigen::Matrix<Scalar, 1, DerivedV::ColsAtCompileTime> q0 = (Scalar(1) - t) * p0 + t * p1;
  Eigen::Matrix<Scalar, 1, DerivedV::ColsAtCompileTime> q1 = (Scalar(1) - t) * p1 + t * p2;
  Eigen::Matrix<Scalar, 1, DerivedV::ColsAtCompileTime> q2 = (Scalar(1) - t) * p2 + t * p3;

  Eigen::Matrix<Scalar, 1, DerivedV::ColsAtCompileTime> r0 = (Scalar(1) - t) * q0 + t * q1;
  Eigen::Matrix<Scalar, 1, DerivedV::ColsAtCompileTime> r1 = (Scalar(1) - t) * q1 + t * q2;

  Eigen::Matrix<Scalar, 1, DerivedV::ColsAtCompileTime> s0 = (Scalar(1) - t) * r0 + t * r1;

  C1.resize(4, C.cols());
  C2.resize(4, C.cols());

  C1.row(0) = p0;
  C1.row(1) = q0;
  C1.row(2) = r0;
  C1.row(3) = s0;

  C2.row(0) = s0;
  C2.row(1) = r1;
  C2.row(2) = q2;
  C2.row(3) = p3;
}

