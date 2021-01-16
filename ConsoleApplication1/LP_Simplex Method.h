#ifndef _DUAL_SIMPLEX_METHOD_H_

#include <iostream>
#include <vector>
#include "Eigen/Dense"

#define IFF 0X7ffffff
using namespace Eigen;

typedef struct MatrixPair {
    MatrixXd Ab;
    MatrixXd An;
}Matrixpair;

double LP_Simplex_Method(std::vector< std::vector<double> > A, std::vector<double> C, std::vector<double>& X, double& value);
Matrixpair get_Mbase(MatrixXd A, VectorXd list, int in, int out, int* & ptr);
int getopt(MatrixXd A, VectorXd b, std::vector<double>& X, double& value);
int Findplus(VectorXd c);
int Findmax(MatrixXd A, int in, VectorXd b);
RowVectorXd getlist(MatrixXd A, VectorXd b);
RowVectorXd randlist(MatrixXd A);
//MatrixXd changeform(std::vector<std::vector<double>> arry, int m, int n);

#endif 