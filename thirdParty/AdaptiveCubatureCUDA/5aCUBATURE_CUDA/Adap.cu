#include "CubaCuda.h"

int main()
{

int dim=1;
int dimf=2560000;
int my_rank=0;
int L=0;

UserInput User;
User.St=new double [dimf];
User.ErrorF=new double [dimf];
User.Nit=10000;
User.TOL=5e-6;
User.TOLr=5e-6;

vector<functor> f(2);


//Call adaptive cubature
for(int i=0;i<2;i++)
{

dim=2-i;

User.xmin=new double [dim];
User.xmax=new double [dim];

for(int j=0;j<dim;j++) 
  {
    User.xmin[j]=0.0;
    User.xmax[j]=100.0;
  }
  
cout<< "Sujo" <<endl;

f[i].function(i);

my_rank=adaptCuba(dim,dimf,f[i],NULL,User,L);

for(int j=0;j<10;j++) {
  cout<< "Converged Integral Result "<< j <<"= "<< setprecision(10) << User.St[j] << endl;}

delete [] User.xmin;
delete [] User.xmax;

sleep(5);

cout<< "Limpo" <<endl;
}



delete [] User.St;
delete [] User.ErrorF;

sleep(5);


  return 0;

}

