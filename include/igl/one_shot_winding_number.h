#ifndef IGL_ONE_SHOT_WINDING_NUMBER
#define IGL_ONE_SHOT_WINDING_NUMBER
#include "igl_inline.h"
#include <Eigen/Core>
#include <vector>
#include <optional>

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

  
template <typename DerivedQ, typename DerivedE, typename DerivedT, typename DerivedS, typename DerivedW>
IGL_INLINE void igl::one_shot_winding_number(
    const Eigen::MatrixBase<DerivedQ> &Q,
    const Eigen::Matrix2<DerivedQ>& E,
    const Eigen::VectorX<DerivedT> & T_sq,
    const Eigen::VectorX<DerivedS> & S,
    const std::optional<Eigen::Matrix2<typename DerivedE::Scalar>>& bounds_opt,
    Eigen::PlainObjectBase<DerivedW> &W);

}

#ifndef IGL_STATIC_LIBRARY
#include "one_shot_winding_number.h"
#endif

#endif

