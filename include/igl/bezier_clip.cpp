#include "bezier_clip.h"
#include "PlainMatrix.h"
#include <cassert>


template<typename DerivedQ, typename DerivedR, typename DerivedC, typename DerivedN>
IGL_INLINE std::pair<std::vector<double>, Eigen::MatrixBase<DerivedN>> bezier_clip(const DerivedQ& origin, const DerivedR& direction, const DerivedC& control_points, double tolerance)
{
  return bezier_clip(origin, direction, control_points, 0.0, 1.0, tolerance);
}


// template<typename DerivedQ, typename DerivedR, typename DerivedC, typename DerivedN>
// IGL_INLINE std::pair<std::vector<double>, Eigen::MatrixBase<DerivedN>> bezier_clip(const DerivedQ& origin, const DerivedR& direction, const DerivedC& control_points, double t_min, double t_max, double tolerance)
// {
//     // Calculate bounding box
//     Matrix2d bounds;
//     bounds << p0.cwiseMin(p1).cwiseMin(p2).cwiseMin(p3),
//               p0.cwiseMax(p1).cwiseMax(p2).cwiseMax(p3);

//     // Check if ray intersects bounding box

//     if (!rayBoxIntersection(origin, direction, bounds)) {
//         return make_pair(VectorXd(), MatrixXd());
//     }
    

//     // If the curve is small enough, consider it an intersection
//     // but also make sure we're far enough from the endpoints points.
//     double tol_sq = tolerance * tolerance;

  
//     if ((p3 - p0).squaredNorm() <= tol_sq) {
//         double t = (tMin + tMax) / 2.0;
       
//         auto [point, tangent] = evaluateBezier(p0, p1, p2, p3, t);

//         Vector2d normal(tangent(1), -tangent(0)); // Perpendicular to the tangent
//         normal.normalize();
//         VectorXd tVec(1);
//         tVec << (origin - point).squaredNorm();
//         MatrixXd normals(1, 2);
//         normals.row(0) = normal;
//         return make_pair(tVec, normals);
//     }

//     // Subdivide the curve
//     double midT = (tMin + tMax) / 2.0;
//     Vector2d p0_left, p1_left, p2_left, p3_left, p0_right, p1_right, p2_right, p3_right;
//     subdivideBezier(p0, p1, p2, p3, midT, p0_left, p1_left, p2_left, p3_left,
//                     p0_right, p1_right, p2_right, p3_right);

//     // Recursively clip the left and right sub-curves
//     auto [tLeft, normalsLeft] = bezier_clip(origin, direction, p0_left, p1_left, p2_left, p3_left, tMin, midT, tolerance);
//     auto [tRight, normalsRight] = bezier_clip(origin, direction, p0_right, p1_right, p2_right, p3_right, midT, tMax, tolerance);

//     // Combine results
//     VectorXd t(tLeft.size() + tRight.size());
//     t << tLeft, tRight;

//     MatrixXd normals(normalsLeft.rows() + normalsRight.rows(), 2);
//     if (normalsLeft.rows() > 0 && normalsRight.rows() > 0) {
//         normals << normalsLeft, normalsRight;
//     } else if (normalsLeft.rows() > 0) {
//         normals = normalsLeft;
//     } else if (normalsRight.rows() > 0) {
//         normals = normalsRight;
//     } else {
//         normals = MatrixXd();
//     }

//     return make_pair(t, normals);
// }



