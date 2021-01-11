#include <iostream>
// #include <math>

using namespace std;
#define MAX_N 200
#define MAX_M 500
int n,m;//size of the problem
double c[MAX_N];
double a[MAX_M][MAX_N];
double b[MAX_M];
double d[MAX_M];
double e[MAX_N];
double delta[MAX_N];
double theta[MAX_N];
int base[MAX_M];//用来记录m个约束，m个基对应n+？个x中的哪些x

double dual_simplex_method();

int main(){
    
    return 0;
}

/*
首先对于给出的矩阵，得到其标准型，即增加若干列，使得不等式变为等式。
之后再将其变成一个有m个初始基的矩阵，即，使得某些列只有一个1，其余全是0.
如此操作请同时给δ（delta）也进行掉（保证基的检验数为0），
且b同步和A矩阵进行了行变换运算。
循环结束的标志：当我判断到所有的b都大于等于0时，说明原始可行，且此时检验数全部小于0，对偶可行。
*/
//入基变量选不重复的。
double dual_simplex_method(){
    double result = 0;
    int base_in_history[MAX_N] = {0};//定义入基的历史个数，入基次数多的尽可能不入基。
    while(1){
        cout<<"检验数:"<<endl;
        for(int i=0;i<n;i++){
            if(delta[i] > 0){
                return -0xffffffff;//对偶问题不可行
            }
            else{
                cout<<delta[i];
            }
        }
        cout<<endl;
        int num = 0;
        for(int i = 0;i<m;i++){
            if(b[i] >= 0){
                num++;
            }
        }
        if(num == m){
            cout<<"原始问题可行"<<endl;
            for(int j=0;j<m;j++){
                result += c[base[j]]*b[base[j]];//当前的x值乘以原来的c值得到结果。
            }
            return result;
        }
        else{//b中一定有小于0的
            int cur_out_base = 0;
            int first_in_base = 0;//这个变量没用，就是为了给smallest_in_base变量做初始化，防止它找不到第一个符合条件的theta。
            int smallest_in_base;
            for(int j=0;j<m;j++){
                if(b[j] < 0&&b[j] < b[cur_out_base]){
                cur_out_base = j;//这里得到的是m个约束中的序号
                }
            }//选定了出基变量。再选入基变量。
            for(int i=0;i<n;i++){
                if(a[cur_out_base][i] < 0){
                    if(first_in_base == 0){
                        smallest_in_base = i;
                        base_in_history[i]++;
                        first_in_base = 1;
                    }
                    theta[i] = delta[i]/a[cur_out_base][i];//不可能会出现除以0的情况，因为考虑的入基，aji都是小于0的。
                    if(theta[i] < theta[smallest_in_base]){
                        //
                        smallest_in_base = i;//更新最小的theta，准备入基。
                        base_in_history[i]++;
                    }
                    else if(theta[i] == theta[smallest_in_base]&&base_in_history[i] < base_in_history[smallest_in_base]){
                        smallest_in_base = i;
                        base_in_history[i]++;
                    }
                }
            }//选定了入基变量，此时可以换基了。
            //换基
            base[cur_out_base] = smallest_in_base;//m个基中的第cur_out_base个基，变成刚刚的入基数。
            double transfer[m];//变换其他m+1个约束，使得刚入基的列满足只有一个1，其他全是0。
            for(int j=0;j<m;j++){
                transfer[j] = a[j][base[cur_out_base]]/a[cur_out_base][base[cur_out_base]];
            }
            for(int j=0;j<m;j++){
                for(int i=0;i<n;i++){
                    if(j == cur_out_base){
                        a[cur_out_base][i] = a[cur_out_base][i]/a[cur_out_base][base[cur_out_base]];//就是当前刚刚入基的列
                    }
                    else{
                        a[j][i] = a[j][i]-a[cur_out_base][i]*transfer[j];
                    }
                }
                if(j == cur_out_base)b[j] = b[j]/a[cur_out_base][base[cur_out_base]];
                else b[j] = b[j]-b[cur_out_base]*transfer[j];
            }//修改所有行列，使得整个矩阵又变成标准形式。
            //修改delta
            for(int i = 0; i < n;i++){
                delta[i] = delta[i] - a[cur_out_base][i]*delta[base[cur_out_base]]/a[cur_out_base][base[cur_out_base]];
            }
        }

    }
    

    
}