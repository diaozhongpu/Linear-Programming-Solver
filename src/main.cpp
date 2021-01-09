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

using namespace std;


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
	for(j=0; j<a.size(); j++)
	{
		cout<<"  ";
		VectorPrint(a[j]);
	}
	cout<<"]";
}


int main(void)
{
	int n, m;
	int i, j;
	double dtmp;
	int itmp;
	vector<double> vdtmp;

	vector<double> c;

	vector< vector<double> > a;
	vector<double> b;
	vector<int> d;
	vector<int> e;

	int k;
	double z;
	vector<double> x;

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
		for(i=0; i<n; i++)
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





	// Call Algorithms





	// Output Display
	switch(k)
	{
		case -1: 
			cout<<"No Solution"<<endl;
			break;
		case 0:
			cout<<"Infinite Solution"<<endl;
			break;
		case 1:
			cout<<"Solution:"<<endl;
			cout<<"Objective Function Optimal: "<<z<<endl;
			cout<<"Variable Values"<<x<<endl;
			break;
		default:
			cout<<"ERROR"<<endl;
	}

	return 0;

}