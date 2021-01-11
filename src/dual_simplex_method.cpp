#include <iostream>
#include <vector>
// #include <math>


using namespace std;




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
double dual_simplex_method(vector< vector<double> > a, vector<double> c, vector<double> &x){
    //c有n个，但是a一行有n+1列。
    int n = c.size();
    int m = a.size();
    vector <double> b;//m个约束的等式右边
    vector <double> delta;//检验数，全部为非正表示对偶可行。
    vector <double> theta;//这个是选择入基变量时的依据，正数theta越小，越优先入基。
    vector <int> base;//用来记录m个基，m个基对应n个x中的哪些x
    double result = 0;
    for(int i =0;i<n;i++){
        delta.push_back(c[i]);
    }//从矩阵中得到delta
    for(int j =0;j<m;j++){
        b.push_back(a[j][n-1]);
    }//从矩阵中得到b
    for(int j =0;j<m;j++){
        base.push_back(n-m+j);
    }//从矩阵中得到base
    cout<<"A ori"<<endl;
    for(int j =0;j<m;j++){
        for(int i =0;i<n;i++){
            cout<< a[j][i] <<" ";
        }
        cout<<endl;
    }
    for(int k=0;k<m;k++){
        vector <double> transfer;
        for(int j=0;j<m;j++){
            transfer.push_back(a[j][base[k]]/a[k][base[k]]);
        }
        for(int j=0;j<m;j++){
            for(int i=0;i<n;i++){
                if(j == k){
                    a[k][i] = a[k][i]/a[k][base[k]];
                }
                else{
                    a[j][i] = a[j][i]-a[k][i]*transfer[j];
                }
            }
            if(j == k)b[j] = b[j]/a[k][base[k]];
            else b[j] = b[j]-b[k]*transfer[j];
        }//修改所有行列
        for(int i = 0; i < n;i++){
            delta[i] = delta[i] - a[k][i]*delta[base[k]]/a[k][base[k]];
        }//修改delta
    }
    cout<<"A"<<endl;
    for(int j =0;j<m;j++){
        for(int i =0;i<n;i++){
            cout<< a[j][i] <<" ";
        }
        cout<<endl;
    }
    vector <int> base_in_history(n,0);//定义入基的历史个数，入基次数多的尽可能不入基。
    while(1){
        //判断对偶问题是否可行
        cout<<"检验数 : "<<endl;
        for(int i=0;i<n;i++){
            if(delta[i] > 0){
                cout<<"dual not works"<<endl;
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
            cout<<"simplex works"<<endl;
            for(int j=0;j<m;j++){
                result += c[base[j]]*b[base[j]];//当前的x值乘以原来的c值得到结果。
            }
            break;
        }
        else{//b中一定有小于0的
            int cur_out_base = 0;
            int first_in_base = 0;//这个变量没用，就是为了给smallest_in_base变量做初始化，防止它找不到第一个符合条件的theta。
            int smallest_in_base = 0;
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
                    theta.push_back(delta[i]/a[cur_out_base][i]);//不可能会出现除以0的情况，因为考虑的入基，aji都是小于0的。
                    if(theta[i] < theta[smallest_in_base]){
                        //
                        smallest_in_base = i;//更新最小的theta，准备入基。
                        base_in_history[i]++;
                        // if(base_in_history[i] >= 100){
                        //     cout<<"too many times,break"<<endl;
                        //     break;
                        // }
                    }
                    else if(theta[i] == theta[smallest_in_base]&&base_in_history[i] < base_in_history[smallest_in_base]){
                        smallest_in_base = i;
                        base_in_history[i]++;
                        // if(base_in_history[i] >= 100){
                        //     cout<<"too many times,break"<<endl;
                        //     break;
                        // }
                    }
                }
                else{// if not taken
                    cout<<"cannot let the simplex works"<<endl;
                    return 0;
                }
            }//选定了入基变量，此时可以换基了。
            //换基
            x[base[cur_out_base]] = 0;//在换基之前，把0赋给x,保证只有基才能影响结果。
            base[cur_out_base] = smallest_in_base;//m个基中的第cur_out_base个基，变成刚刚的入基数。
            double transfer[m];//变换其他m+1个约束，使得刚入基的列满足只有一个1，其他全是0。
            for(int j=0;j<m;j++){
                transfer[j] = a[j][base[cur_out_base]] / a[cur_out_base][base[cur_out_base]];
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
    for(int j =0;j<m;j++){
        a[j][n-1] = b[j];
        x[base[j]] = b[j];
    }//把b赋回去，也把数值交给x。
    return result;
}


