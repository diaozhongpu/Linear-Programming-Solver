#ifndef _DUAL_SIMPLEX_METHOD_H_

#include <iostream>
#include <vector>
#include "Eigen/Dense"

#define IFF 0X7fffffff
using namespace Eigen;


typedef struct MatrixPair {
    MatrixXd Ab;
    MatrixXd An;
}Matrixpair;

double getopt(MatrixXd A, VectorXd b);

#endif 
