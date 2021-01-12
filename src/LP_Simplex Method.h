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

double LP_Simplex_Method(std::vector< std::vector<double> > A, std::vector<double> C);
Matrixpair get_Mbase(MatrixXd A, VectorXd list, int in, int out);
double getopt(MatrixXd A, VectorXd b);
int Findplus(VectorXd c);
int Findmax(MatrixXd A, int in, VectorXd b);
//MatrixXd changeform(std::vector<std::vector<double>> arry, int m, int n);

#endif 