#include <iostream>
#include "Eigen/Dense"
#include <vector>
#include<stdlib.h>
#include<time.h>
#include "LP_Simplex Method.h"

using namespace Eigen;

double LP_Simplex_Method(std::vector< std::vector<double> > A, std::vector<double> C, std::vector<double>& X, double& value)
{
	int m = A.size();
	int	n = A[0].size()-1;
	int ret;
	VectorXd b(m);
	MatrixXd a(m + 1, n);
	int i, j;
	for (i = 0; i < m ; i++) {
		for (j = 0; j < n; j++) {
			a(i+1, j) = double(A[i][j]);
		}
	}
	for (i = 0; i < n; i++) {
		a(0, i) = double(C[i]);
	}
	for (i = 0; i < m; i++) {
		b(i) = double(A[i][n]);
	}
	ret = getopt(a, b, X ,value);
	//for (i = 0; i < a.cols(); i++) {
	//	std::cout << "X[i]" << X[i] << std::endl;
	//}
	return ret;
}

RowVectorXd getlist(MatrixXd A, VectorXd b)
{
	int in = -1, out = -1;							//initialize
	int i = 0, j = 0;
	int* ptr = (int*)malloc(sizeof(int) * A.cols());
	VectorXd b0;
	RowVectorXd Cb, Cn;
	MatrixXd Ab, An, Ab_;
	Matrixpair arr;
	RowVectorXd list = A.row(0);
	srand((unsigned)time(NULL));
	for (i = 0; i < 100; i++) {
		list = randlist(A);
		arr = get_Mbase(A, list, in, out, ptr);	//get the base Matrix a
		Ab = arr.Ab.bottomRows(arr.Ab.rows() - 1);
		An = arr.An.bottomRows(arr.An.rows() - 1);
		Ab_ = Ab.inverse();
		b0 = b;
		b0 = Ab_ * b0;
		std::cout << "in loop to find base for " << i << " times" << std::endl;
		for (j = 0; j < b0.size(); j++) {
			if (!(b0(j) >= 0)) {
				break;
			}
		}
		if (j==b.size()) {
			return list;
		}
		else if((i > (double)A.cols() * (double)A.cols()/1.5) && (i > 10) ){
			list(0) = -1;
			return list;
		}
		else {
			continue;
		}
	}
	list(0) = -1;
	return list;
}

RowVectorXd randlist(MatrixXd A)
{
	int i, j, num;
	RowVectorXd list = A.row(0);
	int n = A.cols();
	for (i = 0; i < A.cols(); i++) {
		list(i) = 0;
	}
	for (j = 0; j < A.rows()-1;) {
		num = rand() % n;
		if (list(num) == 0) {
			list(num) = 1;
			j++;
		}
		else {
			continue;
		}
	}
	return list;
}

int getopt(MatrixXd A, VectorXd b, std::vector<double>& X, double& value)
{
	RowVectorXd Ci = A.topRows(1);
	MatrixXd Ab, An, Ab_;
	Matrixpair arr;
	RowVectorXd c;
	RowVectorXd Cb, Cn;
	int* ptr = (int*)malloc(sizeof(int)*A.cols());
	int ret = 1;		//by default has soluytion
	int i;
	value = 0;		//initialize
	RowVectorXd list = A.row(0);			//vector ,element i !=0 means A(:,i) is part of base Matrix
	if (A.rows()!=(b.size()+1)) {			//check A.rows()==b.size()
		std::cout << "Error0: size doesn't match! " << std::endl;
		exit(0);
	}
	list = getlist(A, b);
	if (list(0) == -1) {
		return -1;
	}
	int in = -1;							//initialize
	int out = -1;
	for (i = 0; i < A.cols(); i++) {
		ptr[i] = i;						//initialize
	}
	while (1) {
		arr = get_Mbase(A, list, in, out, ptr);	//get the base Matrix a
		Ab = arr.Ab.bottomRows(arr.Ab.rows()-1);
		An = arr.An.bottomRows(arr.An.rows() - 1);
		Ab_ = Ab.inverse();
		Cb = arr.Ab.topRows(1);				//get the goal coefficiency
		Cn = arr.An.topRows(1);
		A.block(0, 0, 1, Ab.rows()) = MatrixXd::Zero(1, Ab.rows());
		A.block(0, Ab.rows(), 1, A.cols()-Ab.rows()) = Cn - Cb * Ab_ * An;
		A.block(1, 0, Ab.rows(), Ab.rows()) = MatrixXd::Identity(Ab.rows(), Ab.rows());
		A.block(1, Ab.rows(), Ab.rows(), A.cols() - Ab.rows()) = Ab_ * An;
		//std::cout << "A " << std::endl << A << std::endl;
		//after the calculation, the equation's sequence changes, the base para is always putted to left 
		//std::cout << "A " << std::endl << A << std::endl;
		b = Ab_ * b;
		//std::cout << "b " << std::endl << b << std::endl;
		c = A.row(0);

		for (i = 0; i < A.rows()-1; i++) {
			list(i) = 1;
		}
		for (; i < A.cols(); i++) {
			list(i) = 0;
		}
		double* arry = (double*)malloc(sizeof(double) * A.cols());
		if (arry == NULL) {
			std::cout << "memory allocation failed!" << std::endl;
			exit(5);
		}
		for (i = 0; i < A.rows()-1; i++) {
			arry[ptr[i]] = b(i);
		}
		for (; i < A.cols(); i++) {
			arry[ptr[i]] = 0;
		}
		X.clear();
		value = 0;
		for (i = 0; i < A.cols(); i++) {	
			value += Ci(0, i)* (double)(arry[i]);
			double tmp = (double)(arry[i]);
			X.push_back(tmp);
			//std::cout << "array[i]: " << i << " " << arry[i] << std::endl;
		}
		for (i = 0; i < A.cols(); i++) {
			//std::cout << "X[i]" << X[i] << std::endl;
		}
		//std::cout << "value: " << value << std::endl;
		//std::cout << "c: " << c << std::endl;
		if ((in = Findplus(c))<0 ) {
			//std::cout << "in: " << in << std::endl;
			break;
		}
		if ((out = Findmax(A.bottomRows(A.rows() - 1) , in, b)) < 0) {
			if (out == -1) {
				ret = -1; //no solution
			}
			else {
				ret = 0; //without bound
			}
			break;
		}
		list(in) = 1;
		list(out) = 0;
	}
	//std::cout << "final ret " << ret << std::endl;
	//std::cout << "final value " << value << std::endl;
	return ret;
}

Matrixpair get_Mbase(MatrixXd A, VectorXd list, int in, int out, int*& ptr)
{
	Matrixpair ret;
	int i, j, k, m = A.rows();
	ret.Ab = A.leftCols(A.rows()-1);
	ret.An = A.rightCols(A.cols() - A.rows() + 1);
	
	if ((in >= 0) && (out >= 0)) {								//check if it's valid
		if ((in >= A.cols()) || (out >= A.cols())) {		//check if it's out of range
			std::cout << "Error2: in or out element out of range!" << std::endl;
			exit(2);
		}
		list(in) = 1;
		list(out) = 0;
		//std::cout << "list " << list << std::endl;
	}
	int* tmp = (int*)malloc(sizeof(int) * A.cols());
	for (i = 0; i < A.cols();i++) {
		tmp[i] = ptr[i];
	}
	for (i = 0, j = 0, k = 0; i < list.size(); i++) {
		if (list(i) != 0) {
			ret.Ab.block(0, j, m, 1) = A.block(0, i, m, 1);	//save the col which it's part of base Matrix
			ptr[j] = tmp[i];
			j++;
		}
		else {
			ret.An.block(0, k, m, 1) = A.block(0, i, m, 1);	//save the col which it isn't part of base Matrix
			ptr[k+A.rows()-1] = tmp[i];
			k++;
		}
	}
	if (j != m-1) {			//base size if A.rows()
		std::cout << "Error1: number of base in list goes wrong!" << std::endl;
		exit(1);
	}
	return ret;
}

int Findplus(VectorXd c)
{
	int i, ret=-1;
	for (i = 0; i < c.size(); i++) {
		if (c(i) > 0) {
			ret = i;
			break;
		}
	}
	return ret;
}

int Findmax(MatrixXd A, int in, VectorXd b)
{
	int i, min = IFF, index;
	for (i = 0, index = 0; i < A.rows();i++) {
		if (b(i) < 0) {
			index = -1;
			break;
		}
		if (A(i,in)<=0) {
			continue;
		}
		if (b(i) == 0) {
			continue;
		}
		if ((b(i) / A(i, in)) < min) {
			min = b(i) / A(i, in);
			index = i;
		}
	}
	if (min == IFF) {
		index = -2;
	}
	//std::cout << "index: " << index << std::endl;
	return index;
}