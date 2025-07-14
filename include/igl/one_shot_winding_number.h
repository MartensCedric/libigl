#ifndef IGL_ONE_SHOT_WINDING_NUMBER
#define IGL_ONE_SHOT_WINDING_NUMBER
#include "igl_inline.h"
#include <Eigen/Core>
#include <vector>
namespace igl
{

   template <
    typename DerivedA,
    typename DerivedX,
    typename DerivedW>
  IGL_INLINE void one_shot_winding_number(
    const Eigen::MatrixBase<DerivedA> & A,
    const Eigen::MatrixBase<DerivedX> & X,
    Eigen::PlainObjectBase<DerivedW> & W);

  
  template<typename DerivedQ, typename DerivedC, typename DerivedR, typename DerivedW>
  IGL_INLINE void one_shot_winding_number_cubic_bezier(
    const std::vector<Eigen::MatrixBase<DerivedQ>>& Q,
    const Eigen::MatrixBase<DerivedC>& C,
    Eigen::PlainObjectBase<DerivedW>& W);

}

#ifndef IGL_STATIC_LIBRARY
#include "one_shot_winding_number.h"
#endif

#endif

