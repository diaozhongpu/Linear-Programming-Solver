#pragma once

#include <iostream>
#include <vector>
// #include <math>
using namespace std;




//此函数目标是将一个标准型函数，以及给定的基，进行行变换，得到部分为01矩阵的A。
void row_transformation(vector< vector<double> > a,vector<double> b,vector<double> delta, vector<int> base, int m,int n);


/**
 * 首先对于给出的矩阵，得到其标准型，即增加若干列，使得不等式变为等式。
 * 之后再将其变成一个有m个初始基的矩阵，即，使得某些列只有一个1，其余全是0.
如此操作请同时给δ（delta）也进行掉（保证基的检验数为0），
且b同步和A矩阵进行了行变换运算。
循环结束的标志：当我判断到所有的b都大于等于0时，说明原始可行，且此时检验数全部小于0，对偶可行。
 * 得到的参数：a,c,x,x用于作为输出，同时函数返回值可以直接得到最优解。
 * c相当于检验数，a的最后一列直接赋给b,并且之后一直使用b来进行运算。
 * b最后赋值给x即可。  
 * 入基变量选不重复的。得到解决。目前做法是，记录每个变量入基的次数，优先让入基次数少的入基。
 * 如果还是出现死循环，策略：当他入基次数特别多的时候，输出一下每次的矩阵，以及对应theta，并且使程序停下来。（暂未实现）
*/
double dual_simplex_method(vector< vector<double> > a, vector<double> c, vector<double> &x);
