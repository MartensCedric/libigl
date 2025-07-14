#include "bezier_clip.h"
#include "bezier.h"
#include "ray_box_intersect.h"
#include "PlainMatrix.h"
#include <cassert>


template<typename DerivedQ, typename DerivedR, typename DerivedC, typename DerivedN, typename Scalar>
IGL_INLINE std::pair<std::vector<double>, Eigen::MatrixBase<DerivedN>> igl::bezier_clip(const DerivedQ& origin, const DerivedR& direction, const DerivedC& control_points, Scalar tolerance)
{
  return bezier_clip(origin, direction, control_points, 0.0, 1.0, tolerance);
}

 template<typename DerivedQ, typename DerivedC, typename Scalar>
  std::pair<std::vector<typename DerivedQ::Scalar>, Eigen::MatrixBase<DerivedQ>> igl::bezier_clip(const DerivedQ& origin, const DerivedQ& direction, const DerivedC& control_points, Scalar t_min, Scalar t_max, Scalar tolerance)
  {
    Eigen::Matrix2<typename DerivedQ::Scalar> bounds;

    const auto& p0 = control_points.row(0);
    const auto& p1 = control_points.row(1);
    const auto& p2 = control_points.row(2);
    const auto& p3 = control_points.row(3);

    bounds << p0.cwiseMin(p1).cwiseMin(p2).cwiseMin(p3),
              p0.cwiseMax(p1).cwiseMax(p2).cwiseMax(p3);

 // Check if ray intersects bounding box

    if (!ray_box_intersect(origin, direction, bounds)) {
        return std::make_pair(Eigen::VectorX<Scalar>(), Eigen::MatrixX<DerivedQ>());
    }
    

    // If the curve is small enough, consider it an intersection
    // but also make sure we're far enough from the endpoints points.
    Scalar tol_sq = tolerance * tolerance;

  
    if ((p3 - p0).squaredNorm() <= tol_sq) {
        Scalar t = (t_min + t_max) / 2.0;

        Eigen::Vector2<Scalar> point;
        Eigen::Vector2<Scalar> tangent;
        
        igl::bezier(control_points, t, point, tangent);

        Eigen::Vector2<Scalar> normal(tangent(1), -tangent(0)); // Perpendicular to the tangent
        std::vector<Scalar> t_vec{ (point - origin).norm() };
        DerivedQ normals(1, 2);
        normals.row(0) = normal;
        
        return std::make_pair(t_vec, normals);
      }

    // Subdivide the curve
    Scalar midT = (t_min + t_max) / 2.0;
    DerivedC control_points1;
    DerivedC control_points2;
    subdivideBezier(control_points, midT, control_points1, control_points2);

    // Recursively clip the left and right sub-curves
    auto [tLeft, normalsLeft] = bezier_clip(origin, direction, control_points1, t_min, midT, tolerance);
    auto [tRight, normalsRight] = bezier_clip(origin, direction, control_points2, midT, t_max, tolerance);

    // Combine results
    Eigen::VectorX<Scalar> t(tLeft.size() + tRight.size());
    t << tLeft, tRight;

    Eigen::MatrixX<Scalar> normals(normalsLeft.rows() + normalsRight.rows(), 2);
    if (normalsLeft.rows() > 0 && normalsRight.rows() > 0) {
        normals << normalsLeft, normalsRight;
    } else if (normalsLeft.rows() > 0) {
        normals = normalsLeft;
    } else if (normalsRight.rows() > 0) {
        normals = normalsRight;
    } else {
        normals = Eigen::MatrixX<typename DerivedQ::Scalar>();
    }

    return std::make_pair(t, normals);
  }
