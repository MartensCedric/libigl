#include "one_shot_winding_number.h"
#include "parallel_for.h"
#include "PI.h"
#include <vector>
#include <cassert>
#include <random>

using engine = std::mt19937;

// constexpr double abs_tol = 1e-7;
// constexpr double edge_tol = 1e-7;
// constexpr double martens_tol = edge_tol; 

template <
  typename DerivedA,
  typename DerivedX,
  typename DerivedW>
IGL_INLINE void one_shot_winding_number(
  const Eigen::MatrixBase<DerivedA> & A,
  const Eigen::MatrixBase<DerivedX> & X,  
  Eigen::PlainObjectBase<DerivedW> & W)
{
  assert(A.rows() == X.rows() && "Area and chi matrix must have matching number of rows.");
  assert(A.cols() == X.cols() && "Area and chi matrix must have matching number of columns.");
  W.resize(A.rows());

  for(int i = 0; i < A.rows(); i++)
  {
    W(i) = A.row(i).dot(X.row(i));
  }
}

// VectorXd windingNumberEndpoints(const MatrixXd& C, const MatrixXd& row_query_points, double tolerance) {
//     int wnSize = row_query_points.rows();
//     VectorXd wn = VectorXd::Zero(wnSize);

//     // Extract control points of the Bézier curve
//     Vector2d startPoint = C.row(0);
//     Vector2d endPoint = C.row(3);

//     // Direction vector for ray intersection (horizontal ray)
    
//     Vector2d dir = row_query_points.row(1) - row_query_points.row(0);
//     dir.normalize();
//     Vector2d p0 = C.row(0).transpose();
//     Vector2d p1 = C.row(1).transpose();
//     Vector2d p2 = C.row(2).transpose();
//     Vector2d p3 = C.row(3).transpose();
    
//     Point2D p_start = Point2D{p0(0), p0(1)};
//     Point2D p_end = Point2D{p3(0), p3(1)};
//       Matrix2d bounds;
// 	       bounds << p0.cwiseMin(p1).cwiseMin(p2).cwiseMin(p3),
// 		      p0.cwiseMax(p1).cwiseMax(p2).cwiseMax(p3);

    
//     auto [ts_sq, normals] = rayIntersectBezier(row_query_points.row(0).transpose(), dir,
//                                                 p0, p1, p2, p3, tolerance);
//     std::vector<int> sign(ts_sq.size());
//     for(int k = 0; k < normals.rows(); k++)
//     {
//     	bool same_dir = normals.row(k).dot(dir) > 0.0;
//     	sign[k] = same_dir ? 1 : -1;
//     }
//     Vector3d n(0.0, 0.0, 1.0);
//     Vector3d dir_3d(dir(0), dir(1), 0.0);
//     Vector2d first_q = row_query_points.row(0);  
                                
//     for (int i = 0; i < wnSize; ++i) {
//         Vector2d q = row_query_points.row(i);
        
//         if(!is_inside_boundingbox(bounds, q))
//         {
//            Point2D p{q(0), q(1)};
//            wn(i) = axom::primal::detail::linear_winding_number(p, p_start, p_end, edge_tol);
//         }
//         else
//         {
//                //std::cout << "startPoint: \n" << startPoint << std::endl;
//                //std::cout << "endPoint: \n" << endPoint << std::endl;
//                //std::cout << "q: \n" << q << std::endl;
// 		Vector2d dir_to_start = (startPoint - q).normalized();
// 		Vector2d dir_to_end = (endPoint - q).normalized();
	
// 		double theta_start = acos(dir_to_start.dot(dir));
// 		double theta_end = acos(dir_to_end.dot(dir));
// 		double cosTheta = dir_to_start.dot(dir_to_end);
// 		double thetaRadians = acos(cosTheta);
	 	
// 		Vector3d dir_to_start_3d = Vector3d(dir_to_start(0), dir_to_start(1), 0.0);
// 	 	Vector3d dir_to_end_3d = Vector3d(dir_to_end(0), dir_to_end(1), 0.0);
	 	
// 	 	//std::cout << "dir_to_start: \n" << dir_to_start << std::endl;
// 	 	//std::cout << "dir_to_end: \n" << dir_to_end << std::endl;
// 	 	//std::cout << "cosTheta: " << cosTheta << std::endl;
	 	
// 		//std::cout << "theta_start: " << theta_start << std::endl;
// 		//std::cout << "theta_end: " << theta_end << std::endl;
// 		//std::cout << "thetaRadians: " << thetaRadians << std::endl;
	 	
// 	 	if(dir_to_start_3d.cross(dir_3d).dot(n) > 0.0)
// 	 		theta_start = 2.0 * M_PI - theta_start;

// 		if(dir_to_end_3d.cross(dir_3d).dot(n) > 0.0)
// 	 		theta_end = 2.0 * M_PI  - theta_end; 		
// 	       // Calculate bounding box
	     
// 	       if(dir_to_start_3d.cross(dir_to_end_3d).dot(n) < 0.0)
// 	       {
// 	    		thetaRadians = 2.0 * M_PI  - thetaRadians;
// 	       }


// 	   //std::cout << "theta_start: " << theta_start << std::endl;
// 	   //std::cout << "theta_end: " << theta_end << std::endl;
// 	   //std::cout << "thetaRadians: " << thetaRadians << std::endl;

// 	   double current_t_sq = (q - first_q).squaredNorm();
// 	   int chi = 0;
// 	   for(int k = 0; k < ts_sq.rows(); k++)
// 	   {
// 		if(current_t_sq <= ts_sq(k))
// 		{
// 		 chi += sign[k];
// 		}
// 	   }
	   
// 	   assert(theta_start >= 0.0);
// 	   assert(theta_end >= 0.0);
// 	   assert(thetaRadians >= 0.0);
	   
// 	   int other_chi = chi;
// 	   if(theta_start < theta_end) // inside
// 	   {
// 	    //std::cout << "inside" << std::endl;
// 	    other_chi++;
// 	    //std::cout << "wn(i) = (chi * (2.0 * M_PI - thetaRadians) + other_chi * thetaRadians) / (2.0 * M_PI); " << std::endl;
// 	    //std::cout << chi << " * " << (2.0 * M_PI - thetaRadians) << " + " << other_chi << " * " << thetaRadians << std::endl;
// 	    wn(i) = (chi * (2.0 * M_PI - thetaRadians) + other_chi * thetaRadians) / (2.0 * M_PI); 
	    
// 	   }	   
// 	   else
// 	   {
// 	     //std::cout << "outside" << std::endl;
// 	     other_chi--;
// 	     //std::cout << "wn(i) = ((chi * thetaRadians) + other_chi * (2.0 * M_PI - thetaRadians)) / (2.0 * M_PI); " << std::endl;
// 	     //std::cout << chi << " * " << thetaRadians << " + " << other_chi << " * " << (2.0 * M_PI - thetaRadians) << std::endl;
// 	     wn(i) = ((chi * thetaRadians) + other_chi * (2.0 * M_PI - thetaRadians)) / (2.0 * M_PI); 
	     
// 	   }
//         }	
//     }

//     return wn;
// }
