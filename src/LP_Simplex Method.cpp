/*
author:zzy
说明：这个代码可以独立执行，要求如下：
1.化为等式矩阵
2.变量取址大于等于0
3.求最大值
4.将初始基（对应可行解）传入list（函数增加一个参数）
5.目标不含常数
要求前端转化为相关形式

另：要求输入A，b
A为带目标值的系数矩阵（包含检验数一行）
b为约束等式右边数字
目前简单的例子测试可行
需要测试：
1复杂例子
2无解情况
3无限解情况
4约束有多余（不为行满秩）
*/

#include <iostream>
#include "Eigen/Dense"

#define IFF 0X7fffffff 
using namespace Eigen;


typedef struct MatrixPair {
	MatrixXd Ab;
	MatrixXd An;
}Matrixpair;

Matrixpair get_Mbase(MatrixXd A, VectorXd list, int in, int out);
double getopt(MatrixXd A, VectorXd b);
int Findplus(VectorXd c);
int Findmax(MatrixXd A, int in, VectorXd b);
double getvalue(VectorXd b);

int main()
{
	//example1
	/*MatrixXd A(3,4);
	VectorXd b(2);
	A << 1, 2, 0, 0,
		1, 1, 1, 0,
		0, 1, 0, 1;
	b << 3, 1;*/

	//example2
	/*MatrixXd A(4,5);
	VectorXd b(3);
	A << 0, 1, -2, 0, 0,
		1, -2, 1, 0, 0,
		0, 1, -3, 1, 0,
		0, 1, -1, 0, 1;
	b << 2, 1, 2;*/

	//example3
	MatrixXd A(3,4);
	VectorXd b(2);
	A << -1, 0, -1, 0,
		1, 2, 0, 1,
		0, 1, 2, 0;
	b << 5, 6;

	double result = 0;
	//std::cout << "123" << std::endl;
	result = getopt(A, b);
	std::cout << "The result is: " << result <<std::endl;
	return 0;
}

double getopt(MatrixXd A, VectorXd b)
{
	MatrixXd Ab, An, Ab_;
	MatrixXd tmp;
	Matrixpair arr;
	RowVectorXd c;
	RowVectorXd Cb, Cn;
	double value = 0;
	int i;
	RowVectorXd list = A.row(0);			//vector ,element i !=0 means A(:,i) is part of base Matrix
	//std::cout << A.rows() << " " << b.size() << std::endl;
	if (A.rows()!=(b.size()+1)) {			//check A.rows()==b.size()
		std::cout << "Error0: size doesn't match! " << std::endl;
		exit(0);
	}
	for (i = 0; i < A.cols() - A.rows() + 1; i++) {
		list(i) = 0;					//initialize
	}
	for (; i < A.cols(); i++) {
		list(i) = 1;
	}
	std::cout << "A " << std::endl << A << std::endl;
	//std::cout << "list " << std::endl << list << std::endl;
	//list()
	int in = -1;							//initialize
	int out = -1;
	while (1) {
		arr = get_Mbase(A, list, in, out);	//get the base Matrix a
		//std::cout << "ret.Ab " << std::endl << arr.Ab << std::endl;
		//std::cout << "ret.An " << std::endl << arr.An << std::endl;
		Ab = arr.Ab.bottomRows(arr.Ab.rows()-1);
		An = arr.An.bottomRows(arr.An.rows() - 1);
		//std::cout << "ret.Ab " << std::endl << Ab << std::endl;
		//std::cout << "ret.An " << std::endl << An << std::endl;
		Ab_ = Ab.inverse();
		//std::cout << "Ab " << std::endl << Ab << std::endl;
		//std::cout << "Ab.inverse() " << std::endl << Ab_ << std::endl;
		Cb = arr.Ab.topRows(1);				//get the goal coefficiency
		Cn = arr.An.topRows(1);
		A.block(0, 0, 1, Ab.rows()) = MatrixXd::Zero(1, Ab.rows());
		//std::cout << "A " << std::endl << A << std::endl;
		//std::cout << "Cb.adjoint() " << std::endl << Cb.adjoint()  << std::endl;
		A.block(0, Ab.rows(), 1, A.cols()-Ab.rows()) = Cn - Cb * Ab_ * An;
		//std::cout << "A " << std::endl << A << std::endl;
		A.block(1, 0, Ab.rows(), Ab.rows()) = MatrixXd::Identity(Ab.rows(), Ab.rows());
		//std::cout << "A " << std::endl << A << std::endl;
		A.block(1, Ab.rows(), Ab.rows(), A.cols() - Ab.rows()) = Ab_ * An;
		//after the calculation, the equation's sequence changes, the base para is always putted to left 
		//std::cout << "A " << std::endl << A << std::endl;
		b = Ab_ * b;
		//std::cout << "b " << std::endl << b << std::endl;
		tmp = Cb * b;
		//std::cout << "tmp " << std::endl << tmp << std::endl;
		value += tmp(0);
		//std::cout << "value " << std::endl << value << std::endl;
		c = A.row(0);
		//std::cout << "A " << std::endl << A << std::endl;

		for (i = 0; i < A.rows()-1; i++) {
			list(i) = 1;
		}
		for (; i < A.cols(); i++) {
			list(i) = 0;
		}
		//std::cout << "list " << list << std::endl;

		if ((in = Findplus(c))<0 ) {
			//std::cout << "in: " << in << std::endl;
			break;
		}
		if ((out = Findmax(A.bottomRows(A.rows() - 1) , in, b)) < 0) {
			//std::cout << "out " << out << std::endl;
			break;
		}
		//std::cout << "in: " << in << std::endl;
		//std::cout << "out: " << out << std::endl;
		list(in) = 1;
		list(out) = 0;
		//std::cout << "list " << list << std::endl;
	}
	return value;
}

Matrixpair get_Mbase(MatrixXd A, VectorXd list, int in, int out)
{
	Matrixpair ret;
	int i, j, k, m = A.rows();
	ret.Ab = A.leftCols(A.rows()-1);
	ret.An = A.rightCols(A.cols() - A.rows() + 1);
	//std::cout << "list " << std::endl << list << std::endl;
	//std::cout << "ret.Ab " << std::endl << ret.Ab << std::endl;
	//std::cout << "ret.An " << std::endl << ret.An << std::endl;
	if ((in >= 0) && (out >= 0)) {								//check if it's valid
		if ((in >= A.cols()) || (out >= A.cols())) {		//check if it's out of range
			std::cout << "Error2: in or out element out of range!" << std::endl;
			exit(2);
		}
		list(in) = 1;
		list(out) = 0;
		std::cout << "list " << list << std::endl;
	}
	
	for (i = 0, j = 0, k = 0; i < list.size(); i++) {
		if (list(i) != 0) {
			//std::cout << list << std::endl;
			ret.Ab.block(0, j, m, 1) = A.block(0, i, m, 1);	//save the col which it's part of base Matrix
			j++;
			//std::cout << "ret.Ab " << std::endl << ret.Ab << std::endl;
			//std::cout << "i " << std::endl << i << std::endl;
		}
		else {
			ret.An.block(0, k, m, 1) = A.block(0, i, m, 1);	//save the col which it isn't part of base Matrix
			//std::cout << "ret.An " << std::endl << ret.An << std::endl;
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
			//return -1;		//no solution
			std::cout << "Error3: no solution!" << std::endl;
			exit(3);
		}
		if (A(i,in)<=0) {
			//return -2;		//without bound
			
			continue;
		}
		if (b(i) == 0) {
			continue;
		}
		//std::cout << "(b(i) / A(i, in): " << b(i) / A(i, in) << std::endl;
		if ((b(i) / A(i, in)) < min) {
			min = b(i) / A(i, in);
			index = i;
			//std::cout << "i: " << i << std::endl;
		}
	}
	if (min == IFF) {
		std::cout << "Error4: without bound!" << std::endl;
		exit(4);
	}
	//std::cout << "hello" << std::endl;
	return index;
}