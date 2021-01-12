#include <iostream>
#include "Eigen/Dense"
 
using namespace Eigen;
 
int main()
{
  int Dynamic = 2;
  Matrix2d a;
  //Coefficient accessors
  a(0,0) = 3;
  a(1,0) = 2.5;
  a(0,1) = -1;
  a(1,1) = a(1,0) + a(0,1);
  MatrixXd b(Dynamic,Dynamic);
  MatrixXd bi(2,2);
  //Comma-initialization
  b << 2, 3,
       1, 4;
  std::cout << "a + b =\n" << a + b << std::endl;
  std::cout << "a - b =\n" << a - b << std::endl;
  std::cout << "a_11 =\n" << a(1,1) << std::endl;
  std::cout << "Doing a += b;" << std::endl;
  a += b;
  std::cout << "Now a =\n" << a << std::endl;
  Vector3d v(1,2,3);
  Vector3d w(1,0,0);
  std::cout << "-v + w - v =\n" << -v + w - v << std::endl;

  bi=b.inverse();
  std::cout << "NowNow bi =\n" << bi << std::endl;

  std::cout << "NowNow I =\n" << b*bi << std::endl;

}