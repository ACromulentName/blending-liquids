#define EIGEN_USE_MKL_ALL
#include <Eigen/Sparse>
#include <Eigen/Dense>

#define SINGLE_PRECISION

#ifdef SINGLE_PRECISION
typedef float Float;
typedef Eigen::VectorXf VectorX;
typedef Eigen::MatrixXf MatrixX;
typedef Eigen::Matrix3f Matrix3;
typedef Eigen::Matrix4f Matrix4;
#else
typedef double Float;
typedef Eigen::VectorXd VectorX;
typedef Eigen::MatrixXd MatrixX;
typedef Eigen::Matrix3d Matrix3;
typedef Eigen::Matrix4d Matrix4;
#endif

typedef Eigen::Matrix<Float, 4, 1> Vector4;
typedef Eigen::Matrix<Float, 3, 1> Vector3;
typedef Eigen::Matrix<Float, 2, 1> Vector2;

namespace Globals
{
	extern Float AvgEdgeLength;
	extern const int InvalidNeighbour;	
}

extern Eigen::IOFormat CommaFmt; 

#define ENABLE_PARALLEL 1