#ifndef CubaCuda_H
#define CubaCuda_H


#include <iterator>
#include <string>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <vector>


typedef double (*func) (const double *,int, int);

#define HDF __host__ __device__
#define DF __device__


HDF double F(const double *y, int index,int id); 


class functor
{
public:
  void function(int index) { index_ = index ;}
   HDF
    double operator()(const double *y,int index,int id)
    {
        return F(y,index_,id);
    } 

 private:
  int index_;  
        
};

using namespace std;



#define dimk(dim) (dim*(2*dim*(dim-1)))
#define diml(dim) (4*dim*dim)
#define dimm(dim) ((1U << (dim))*dim)



const  double alpha=0.3585685828003180919906451539079374954541;
const  double beta=0.9486832980505137995996680633298155601160;
const  double gammas=0.6882472016116852977216287342936235251269;


unsigned ls0(unsigned n);

void evalR_Rfs(int dim,double *pts,  double *p, const double *c, const double *r);

void evalRR0_0fs(int dim,double *pts, double *p, const double *c, const double *r);

void evalR0_0fs4d(int dim,double *pts, double *p, const double *c,
			 const double *r1, const double *r2);


void points(int dim,double *pts1,double *pts2,double *pts3);



class Data
{   
 public:
 
  //---------------PINNED MEMORY-------------//
   double *x0;
   double *x1;
   double *x2;
   double *x3;
   double *diff;
   double *ERROR;
   double *S7;
//---------------PINNED MEMORY-------------//  
 
}; 


class DataDU
{   
 public:
 
  double *d_ERROR;
  double *d_S7;
 
 
};


class DataD: public DataDU
{   
 public:
 
  //---------------DATA TO DEVICE-------------//
   double *d_S5;
   double *d_x1;
   double *d_x2;
   double *d_x3;
   double *d_diff;
   double *d_W1;
   double *d_W;
   double *d_x0;
};


class UserInput
{   
 public:
 
  //---------------User Input-------------//
   int Nit;
   double TOL;
   double TOLr;
   double *xmin;
   double *xmax; 
   double *St;
   double *ErrorF;   
}; 




class Region
{   
 public:
// Input Geometric coordinate;
   vector<double> a;
  //  double a[8];
   vector<double> b;
  // double b[8];
// OutPut Information;
// Integral Result
   vector<double> S;
   vector<double> Errorabs;

    int split;
    double ErrorMax;
// constants
  void  Declare(int dim, int dimf)
   {
    S.resize(dimf);
    Errorabs.resize(dimf);
    a.resize(dim);
    b.resize(dim);   
   }
   
// Domain Definition
   void Regions(int dim, double *dx, double *d);   
// Cubature Member fuctions
   void CubaCuda(int dim, int dimf,functor f, void * par,double *pts1,double *pts2,double *pts3, double *W, double *W1,int my_rank,int proc, const Data & dat);
   void Cubature(int dim,int dimf, double *dx, double *d,double *pts1,double *pts2,double *pts3, functor f,void * par,double *S5,double *W, double *W1,int my_rank,int proc, const Data& dat);
//Quadrature Member functions
   void Quadrature(int dim,int dimf, double *dx, double *d, functor f,void * par,int my_rank,int proc, const Data& dat);  

  
};

int adaptCuba(int dim,int dimf, functor f, void * par, const UserInput& User, int& L);


#endif 
