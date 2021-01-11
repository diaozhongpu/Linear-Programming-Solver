/*
对偶单纯形算法
输入格式：
A：单纯形表矩阵，各变量的系数和常数 M * N+1 维
a_11 a_12 ··· a_1N b_1
a_21 a_22 ··· a_2N b_2
 ·     ·        ·   ·
 ·     ·        ·   ·
 ·     ·        ·   ·
a_M1 a_M2 ··· a_MN b_M

C：目标函数系数向量 N 维
Lambda：检验数向量 N 维
base_V：基向量下标 M 维
B：基向量在目标函数中对应的系数 M 维
goal：目标值

flygeda 2017.12.16
*/
#include <stdio.h>
#include <iostream>

#define N 6  
#define M 3 //N个变量，M个约束条件

using namespace std;

double A[M][N + 1] = { { 3, -1, 1, -2, 0, 0, -3},
                    { 2, 1, 0, 1, 1, 0, 4 },
                    { -1, 3, 0, -3, 0, 1, 12 } };
double C[N] = { -7, 7, -2, -1, -6, 0 };//目标函数系数
//double Goal[N] = { 0, 0, 0, 0, 0, 0};//解向量
double Lambda[N] = { 0, 0, 0, 0, 0, 0 };//检验数
int base_V[M] = { 2, 4, 5 };//基向量，初始基向量为 矩阵A 的第2，4，5列
double B[M] = { 0, 0, 0 };
double goal;//目标值

bool is_base(int v)
{
    for (int i = 0; i < M; i++)
        if (base_V[i] == v)
            return true;
    return false;
}

double dual()
{
    // 给出矩阵A的形式，便于分析
    for (int i = 0; i < M; i++) 
    {
        for (int j = 0; j < N + 1; j++)
            cout << " " << A[i][j];
        cout << endl;
    }

    while(1)
    {
        for (int i = 0; i < M; i++)
            B[i] = C[base_V[i]];
        //计算检验数
        cout << "检验数C:" << endl;
        for (int i = 0; i < N; i++)
        {
            double lam_temp = 0;
            for (int j = 0; j < M; j++)
                lam_temp += B[j] * A[j][i];
            Lambda[i] = C[i] - lam_temp;

            if (Lambda[i] < 0)
            {
                cout << "can not be soluted by dual method" << endl;
                return 0;
            }
            cout << Lambda[i] << ' ';
        }
        cout << endl;
        //计算目标值
        double goal_temp = 0;
        for (int i = 0; i < M; i++)
            goal_temp += B[i] * A[i][N];
        goal = goal_temp;
        cout << "目标值" << goal << endl;

        int temp = 0;
        double min_b = 1000;
        int out_base = 0;
        int row = 0;
        cout << "常数列" << endl;
        for (int i = 0; i < M; i++)
        {
            cout << A[i][N] << ' ';
            if (A[i][N] < 0)//判断常数列
                temp = 1;
            if (A[i][N] < min_b)
            {
                min_b = A[i][N];
                out_base = base_V[i]; //确定出基变量下标
                row = i;
            }
        }
        cout << endl;
        //取得最优值条件
        if (temp == 0)
            break;
        //未达到最优条件
        double min_abs = 1000;
        int in_base = 0;
        for (int i = 0; i < N; i++)
        {

            if (!is_base(i) && A[row][i] < 0)
            {
                double temp_min = abs(Lambda[i] / A[row][i]);
                if (temp_min < min_abs)
                {
                    min_abs = temp_min;//确定入基变量
                    in_base = i;
                }
            }
        }
        base_V[row] = in_base;
        //矩阵行变换，换基运算
        for (int i = 0; i < M; i++)
        {
            temp = A[row][in_base];
            for (int j = 0; j < N + 1; j++)
            {
                if (i != row)
                    A[i][j] += -A[row][j] / temp;
                else
                    A[i][j] /= temp;
            }
        }
        //cout << "in_base:" << in_base << "\trow:" << row << endl;
        cout << endl;
    }
    return goal;
}

int main()
{
    double g = dual();
    cout << endl;
    cout << "最优值" << g << endl;
    //cout << A[0][6];
    return 1;
}