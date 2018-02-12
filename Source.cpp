#include <math.h>
#include <iostream>
#include <iomanip>
#include "Forsythe.h"

using namespace std;

double** build(double **A, double EPS);
double norma(double **A);
double** function(double **A, double **R);
double** obrat(double **A);
double** multiplier(double **A, double **B);

const int size = 5;

void main(){
	/*              test             */
	/*double *test1 = new double [size];
	double *test2 = new double [size];


	for(int i=0; i<size; i++){test1[i]=i;}
	cout<<"test1:"<<endl;
	for(int i=0; i<size; i++){cout<<setprecision(5)<<setw(10)<<test1[i]<<" ";}
	cout<<endl<<endl<<"test2=test1 :"<<endl;
	test2=test1;
	for(int i=0; i<size; i++){cout<<setprecision(5)<<setw(10)<<test2[i]<<" ";}
	cout<<endl<<endl<<"for(int i=0; i<size; i++){test2[i]=test1[i];} :"<<endl;
	for(int i=0; i<size; i++){test2[i]=test1[i];}
	for(int i=0; i<size; i++){cout<<setprecision(5)<<setw(10)<<test2[i]<<" ";}


	double **test3 = new double*[size];
	double **test4 = new double*[size];
	for (int i=0;i<size;i++)
	{
		test3[i]=new double [size];
		test4[i]=new double [size];
	} 

	cout<<endl<<endl<<"test3:"<<endl;
	for(int i = 0; i < size; i++){
	 for(int j = 0; j < size; j++){
		test3[i][j] = i;
		cout<<test3[i][j]<<" ";
	 }
	 cout<<endl;
	}

	
	cout<<endl<<"test4=test3 :"<<endl;
	test4=test3;
	for(int i = 0; i < size; i++){
	 for(int j = 0; j < size; j++){
		cout<<test4[i][j]<<" ";
	 }
	 cout<<endl;
	}
*/
	/*       end of test             */

	cout.fixed;

	double **A = new double*[size];
	double **R = new double*[size];
	for (int i=0;i<size;i++)
	{
		R[i]=new double [size];
		A[i]=new double [size];
	} 
	
    double norm, norm1, norm2,  EPS, cond;
    EPS=0.001;

	for(int i=0; i<3; i++){
		//cout<<"                               Matrix A:"<<endl;
		A=build(A,EPS);
		R=function(A,R);
		norm=norma(R);
		//cout.scientific;
		cout<<endl<<"Norm: "<<norm<<endl<<"EPS: "<<setprecision(8)<<EPS<<endl;
		EPS=EPS/10;
	}
	//**********************************//
	EPS=0.000001;
	//cout<<"                               Matrix A:"<<endl;
	build(A,EPS);
    norm1=norma(R);
	obrat(A);
	norm2=norma(R);
	//cond_MULT=norm1*norm2;
	cout<<endl<<"EPS: "<<setprecision(8)<<EPS<<endl<<"Norm: "<<norm1<<" | Obrat norm: "<<norm2<<" | cond: "<<cond<<endl;

}

//Построение матрицы
double** build(double **A, double EPS){
      double tmp;

	  for(int i = 0; i<size; i++){
		  for(int j = 0; j<size; j++){
			tmp=j+1;
            A[i][j]=(1.0+cos(tmp))/(sin(tmp)*sin(tmp));
			if(j==4) {A[i][j]=(1.0+cos(1.0))/((sin(1.0+EPS)*sin(1.0+EPS)));}
			double power = A[i][j];
            for(int involution = 0; involution<i; involution++){
				A[i][j]=A[i][j]*power;
			}

			//cout<<setprecision(5)<<setw(10)<<A[i][j]<<" |  ";
		  }
		  //cout<<endl;
	  }
	  return A;
}

//нахождение нормы
double norma(double **A){
     double norm=0;
	 double tmp=0;
 
	  for(int i = 0; i<size; i++){
		  tmp=0;
		  for(int j = 0; j<size; j++){
			  tmp=tmp+abs(A[j][i]);
		  }
		 if(tmp>norm) {norm=tmp;}
      }
	  return norm;
}

double** function(double **A, double **R){
	double** result= new double*[size];
	for (int i=0;i<size;i++)
	{
		result[i]=new double [size];
	} 

	for(int i = 0; i < size; i++){
	 for(int j = 0; j < size; j++){
		result[i][j] = A[i][j];
	 }
	}

	result=obrat(result);

	//cout<<endl<<endl<<"E?!"<<endl;
	////cout.fixed;
	//R = multiplier(result, A);
	////cout.scientific;

	//multiplier(A, result);
	for(int i = 0; i<size; i++){
		R[i][i]--;
	}
	cout<<"                               Matrix R"<<endl;
	for(int i=0; i<size; i++){
		for(int j=0; j<size; j++){
			cout<<setprecision(5)<<setw(10)<<R[i][j]<<" | ";
		}
		cout<<endl;
	}

	return R;
} 

//Нахождение обратной с помощью decomp и solve
double** obrat(double **A){
	double** result= new double*[size];
	for (int i=0;i<size;i++)
	{
		result[i]=new double [size];
	} 
	double cond;
	int ipvt[size];
	int k=0;
	double newA[size*size];
	for(int i=0; i<size; i++){
		for(int j=0; j<size; j++){
			newA[k]=A[i][j];
			k++;
		}
	}
	for(int i=0; i<size*size; i++){if(i%5==0){cout<<endl;}cout<<newA[i]<<" ";}
	Decomp(size, newA, &cond, ipvt);
	cout<<endl<<endl;
	for(int i=0; i<size*size; i++){if(i%5==0){cout<<endl;}cout<<newA[i]<<" ";}
	double B[size];
	
	for(int i=0; i<size; i++){
		for(int k=0; k<size; k++){B[k]=0;}
		B[i]=1;
		Solve(size, newA, B, ipvt);
		for(int j=0; j<size; j++){
			result[j][i]=B[j];
		}
	}
	return result;
}

//Умножение матриц
double** multiplier(double **A, double **B){
	double** result= new double*[size];
	for (int i=0;i<size;i++){
		result[i]=new double [size];
	}

	for(int i=0; i<size; i++){
		for(int j=0; j<size; j++){
			result[i][j]=0;
		}
	}

	for(int i=0; i<size; i++){
		for(int j=0; j<size; j++){
			for(int k=0; k<size; k++){
				result[i][j]+=A[i][k]*B[k][j];
			}
		}
	}

	//cout<<"                               Multiplier: "<<endl;
	//for(int i=0; i<size; i++){
	//	for(int j=0; j<size; j++){
	//		//if(result[i][j]<0.0000000001){cout<<setprecision(5)<<setw(10)<<"0 | ";}
	//		//else {cout<<setprecision(5)<<setw(10)<<result[i][j]<<" | ";}
	//		//cout<<setprecision(5)<<setw(10)<<result[i][j]<<" | ";
	//	}
	//	//cout<<endl;
	//}
	return result;
}

