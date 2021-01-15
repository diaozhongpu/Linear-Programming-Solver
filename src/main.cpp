/*
 * Linear Programming Solver Using Simplex Method and Dual Simplex Method
 * By Yanbo Jiang, Zhongpu Diao, and Ziyue Zhan
 *
 * Frontend by Zhongpu Diao
 *
 */

#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "dual_simplex_method.h"
#include "LP_Simplex Method.h"

using namespace std;

void VectorPrintI(vector<int> path)
{
	std::vector<int>::const_iterator i;
	cout<<"[ ";
	for (i = path.begin(); i != path.end(); ++i)
    	cout << *i << ' ';
    cout<<"]"<<endl;
}

void VectorPrint(vector<double> path)
{
	std::vector<double>::const_iterator i;
	cout<<"[ ";
	for (i = path.begin(); i != path.end(); ++i)
    	cout << *i << ' ';
    cout<<"]"<<endl;
}

void MatrixPrint(vector< vector<double> > a)
{
	int j;
	cout<<"["<<endl;
	for(j=0; j < (int)a.size(); j++)
	{
		cout<<"  ";
		VectorPrint(a[j]);
	}
	cout<<"]"<<endl;
}


int main(void)
{
	int n, m;
	int i, j, iXInR, t;
	double dtmp;
	int itmp;
	vector<double> vdtmp;

	vector<double> c;

	vector< vector<double> > a;
	vector<double> b;
	vector<int> d;
	vector<int> e;

	int k = 0;
    int tempk = 0;
	vector<double> realx;
	vector<double> finx;
	double finopt=INFINITY; // min or max?
	vector<int> finE;
	vector<double> x;
	double opt;
    int useSimplex; // 0 for dual, 1 for simplex
    cout<<"input 0 for dual, 1 for simplex"<<endl;
    cin>>useSimplex;
    useSimplex=useSimplex?1:0;
    
    
    
	cout<<"Please input"<<endl;

	// Input Process
	cin>>n;
	cin>>m;
	
	for(i=0; i<n; i++)
	{
		cin>>dtmp;
		c.push_back(dtmp);
	}


	for(j=0; j<m; j++)
	{
		vdtmp.clear();
		for(i=0; i<n; i++) // a
		{
			cin>>dtmp;
			vdtmp.push_back(dtmp);
		}
		a.push_back(vdtmp);

		cin>>dtmp;
		b.push_back(dtmp);

		cin>>itmp;
		d.push_back(itmp);

	}

	for(i=0; i<n; i++)
	{
		cin>>itmp;
		e.push_back(itmp);
	}

    cout<<"a:"; MatrixPrint(a);
    cout<<"b:"; VectorPrint(b);
    cout<<"c:"; VectorPrint(c);
    cout<<"d:"; VectorPrintI(d);
    cout<<"e:"; VectorPrintI(e);



	// 丝考虑超定问题
	// 将等弝放到最剝面
	int artVarNum=0; // number of artificial variables
	vector<int>::iterator dnth;
	vector< vector<double> >::iterator anth;
    vector<double>::iterator bnth;
	vector<int> baseMap;
	for(j=0; j<(int)d.size(); j++)
	{
		if(d[j]==1||d[j]==-1)
		{
			artVarNum++;
		}
		else
		{
            if(useSimplex==0){
            
                baseMap.insert(baseMap.begin(), j);
                d.insert(d.begin(), e[j]);
                dnth = d.begin() + j+1;
                d.erase(dnth);
                a.insert(a.begin(), a[j]);
                anth = a.begin() + j+1;
                a.erase(anth);
                b.insert(b.begin(), b[j]);
                bnth = b.begin() + j+1;
                b.erase(bnth);
            }
            else
            {
            
                d[j]=1;
                artVarNum++;
                a.push_back(a[j]);
                b.push_back(b[j]);
                d.push_back(-1);
                m++;
             
            }
            //artVarNum++;
		}
	}

	//TODO: choose the right base



	// standardize 

	// divide x \in R: e 1, >=0; 0, R; -1, <=0
	vector<int> xInR;
	for(i=0; i<e.size(); i++)
	{
		if(e[i]==0)
		{
			xInR.push_back(i);
		}
	}
    cout<<"xInR:"; VectorPrintI(xInR);
	vector< vector<double> > A;
	vector<double> C;
	vector<int> E;
	if((int)xInR.size()==0)
	{
		// change x in R into 2 seperate x>=0 and x<=0
		// change e
		E=e;//copy


		// add artificial variables +-a.r.; not for == 
		for(i=0; i<c.size(); i++)
		{
			C.push_back(-c[i]);
		}
		//C=c;

		for(j=0; j<m-artVarNum; j++)
		{
			vdtmp.clear();
			vdtmp=a[j];
			for(i=0; i<artVarNum; i++)
			{
				vdtmp.push_back(0);
			}
            vdtmp.push_back(b[j]);
			A.push_back(vdtmp);
		}

		for(j=m-artVarNum; j<m; j++)
		{
			vdtmp.clear();
			vdtmp=a[j];
			for(i=0; i<artVarNum; i++) // a
			{
				if(i==(j-(m-artVarNum)))
				{
					vdtmp.push_back(1);
				}
				else
				{
					vdtmp.push_back(0);
				}
			}
            vdtmp.push_back(b[j]);
			A.push_back(vdtmp);

			C.push_back(0);
            E.push_back(-d[j]); //!!
		}

		// 行坘�?==为基
        /*
        int nE=(int)E.size();
		for(j=1; j<=m-artVarNum; j++)//m-artVarNum-1
		{
			for(i=0; i<(int)A[j].size(); i++)
			{
				if(A[m-artVarNum-j][nE-artVarNum-j]==0)
				{
					printf("ZERO PIVOT!");
                    int fixed=0;
                    // try to find the right base
                    for(t=0; t<m-artVarNum; t++)
                    {
                        
                        if(A[t][nE-artVarNum-j]!=0)
                        {
                            //swap t and m-artVarNum-j
                            if(t<m-artVarNum-j)
                            {
                                anth=A.begin()+t;
                                A.insert(anth, A[m-artVarNum-j]);
                                anth=A.begin()+m-artVarNum-j+1;
                                A.insert(anth, A[t+1]);
                                anth=A.begin()+m-artVarNum-j+2;
                                A.erase(anth);
                                anth=A.begin()+t+1;
                                A.erase(anth);
                                fixed=1;//fixed
                                break;
                            }
                            
                        }
                    }
                    if(fixed)
                    {
                        printf("FIND NONZERO PIVOT!");
                        break;
                    }
                    if(t>=m-artVarNum)
                    {
                        printf("FAILED TO FIND ZERO PIVOT!");
                        exit(1);
                    }
				}
				A[m-artVarNum-j][i]/=A[m-artVarNum-j][nE-artVarNum-j];
			}
			for(i=0; i<m; i++)
			{
				if((m-i)!=j)
				{
					for(t=0; t<(int)A[i].size(); t++)
					{
						A[i][t]=A[i][t]-A[m-artVarNum-j][t]*A[i][nE-artVarNum-j];//A[j][j]==1
					}
				}
			}
		}
         */

		// change variable sign
		for(i=0; i<(int)E.size(); i++)//n+artVarNum
		{
			if(E[i]==-1)
			{
				for(j=0; j<m; j++)
				{
					A[j][i]=-A[j][i];
				}
				C[i]=-C[i];
			}
		}



		// Call Algorithms


		/* M constraints, N variables (with artificial, etc)
		 * A (M)*(N+1) matrix, first N column is variable, of which last M is selected as base variable
		 * C (N) objective function co
		 *
		 */
        for(i=0; i<n+artVarNum; i++)
        {
            finx.push_back(-1);
        }
        printf("Single");
        cout<<"A:"; MatrixPrint(A);
        cout<<"C:"; VectorPrint(C);
		k=dual_simplex_method(A, C, finx, finopt);

		finE=E;


	}
	else
	{
		for(iXInR=0; iXInR<(2<<((int)xInR.size())/2); iXInR++)
		{
			// change x in R into 2 seperate x>=0 and x<=0
			// change e
			E=e;//copy
			if((int)xInR.size()>=1)
			{
				for(j=0; j<(int)xInR.size(); j++)
				{
					if((2*iXInR/(2<<j))%2==1)
					{
						//printf("j=%d, iXInR=%d, (2<<j)/2=%d, E[%d]=1\n", j, iXInR, (2<<j)/2, xInR[j]);
						E[xInR[j]]=1;
					}
					else
					{
						//printf("j=%d, iXInR=%d, (2<<j)/2=%d, E[%d]=-1\n", j, iXInR, (2<<j)/2, xInR[j]);
						E[xInR[j]]=-1;
					}
				}
			}



			// add artificial variables +-a.r.; not for == 

            for(i=0; i<c.size(); i++)
            {
                C.push_back(-c[i]);
            }

            A.clear();
			for(j=0; j<m-artVarNum; j++)
			{
				vdtmp.clear();
				vdtmp=a[j];
				for(i=0; i<artVarNum; i++)
				{
					vdtmp.push_back(0);
				}
                vdtmp.push_back(b[j]);
				A.push_back(vdtmp);
			}

			for(j=m-artVarNum; j<m; j++)
			{
				vdtmp.clear();
				vdtmp=a[j];
				for(i=0; i<artVarNum; i++) // a
				{
					if(i==(j-(m-artVarNum)))
					{
						vdtmp.push_back(1);
					}
					else
					{
						vdtmp.push_back(0);
					}
				}
                vdtmp.push_back(b[j]);
				A.push_back(vdtmp);

				C.push_back(0);
				E.push_back(-d[j]); //!!
			}

			// 行坘�?==为基
            /*
            int nE=(int)E.size();
            for(j=1; j<=m-artVarNum; j++)//m-artVarNum-1
            {
                for(i=0; i<(int)A[j].size(); i++)
                {
                    if(A[m-artVarNum-j][nE-artVarNum-j]==0)
                    {
                        printf("ZERO PIVOT!");
                        int fixed=0;
                        // try to find the right base
                        for(t=0; t<m-artVarNum; t++)
                        {
                            
                            if(A[t][nE-artVarNum-j]!=0)
                            {
                                //swap t and m-artVarNum-j
                                if(t<m-artVarNum-j)
                                {
                                    anth=A.begin()+t;
                                    A.insert(anth, A[m-artVarNum-j]);
                                    anth=A.begin()+m-artVarNum-j+1;
                                    A.insert(anth, A[t+1]);
                                    anth=A.begin()+m-artVarNum-j+2;
                                    A.erase(anth);
                                    anth=A.begin()+t+1;
                                    A.erase(anth);
                                    fixed=1;//fixed
                                    break;
                                }
                                
                            }
                        }
                        if(fixed)
                        {
                            printf("FIND NONZERO PIVOT!");
                            break;
                        }
                        if(t>=m-artVarNum)
                        {
                            printf("FAILED TO FIND ZERO PIVOT!");
                            exit(1);
                        }
                    }
                    A[m-artVarNum-j][i]/=A[m-artVarNum-j][nE-artVarNum-j];
                }
                for(i=0; i<m; i++)
                {
                    if((m-i)!=j)
                    {
                        for(t=0; t<(int)A[i].size(); t++)
                        {
                            A[i][t]=A[i][t]-A[m-artVarNum-j][t]*A[i][nE-artVarNum-j];//A[j][j]==1
                        }
                    }
                }
            }
             */

			// change variable sign
			for(i=0; i<(int)E.size(); i++)//n+artVarNum
			{
				if(E[i]==-1)
				{
					for(j=0; j<m; j++)
					{
						A[j][i]=-A[j][i];
					}
					C[i]=-C[i];
				}
			}



			// Call Algorithms


			/* M constraints, N variables (with artificial, etc)
			 * A (M)*(N+1) matrix, first N column is variable, of which last M is selected as base variable
			 * C (N) objective function co
			 *
			 */
			
            // k
            for(i=0; i<n+artVarNum; i++)
            {
                x.push_back(-1);
            }
            printf("Multi");
            cout<<"A:"; MatrixPrint(A);
            cout<<"C:"; VectorPrint(C);
            if(finopt==INFINITY)
            {
                k=dual_simplex_method(A, C, x, opt);
                finopt=opt;
                finx=x;
                finE=E;
            }
            else
            {
                tempk=dual_simplex_method(A, C, x, opt);
                if(tempk==1)
                {
                    k=1;
                    
                    if(opt<finopt)//update optimal if has x in R
                    {
                        finopt=opt;
                        finx=x;
                        finE=E;
                    }
                }
            }
            
			
			


		}
	}
	

	// Output Display
	// merge opt
	// fix variable sign
	vector<double>::iterator realxnth;
	switch(k)//TODO:k
	{
		case -1: 
			cout<<"No Solution"<<endl;
			break;
		case 0:
			cout<<"Infinite Solution"<<endl;
			break;
		case 1:
			//for(i=n-artVarNum; i<n; i++)
            for(i=0; i<n-(m-artVarNum); i++)
			{
				realx.push_back(finx[i]*finE[i]);
			}
			for(i=0; i<n; i++)
			{
				for(j=0; j<m-artVarNum; j++)
				{
					if(baseMap[j]==i)
					{
						realxnth=realx.begin()+i;
						realx.insert(realxnth, finx[j]*finE[j]);
					}
				}
			}
			/*
            t=n-1;
			for(i=0; i<artVarNum; i++)
			{
				while(t>=0&&d[t]!=0)//n=d.size()
				{
					t--;
				}
				if(t>=0)
				{
					realxnth=realx.begin()+t;
					realx.insert(realxnth, finx[i]*finE[i]);
					t--;
				}
			}
			*/
			//i==artVarNum, t==-1/last d[i]
			cout<<"Solution:"<<endl;
			cout<<"Objective Function Optimal: "<<-finopt<<endl;
			cout<<"Variable Values: "<<endl;
			VectorPrint(realx);
			cout<<endl;
			break;
		default:
			cout<<"ERROR"<<endl;
	}

	return 0;

}
