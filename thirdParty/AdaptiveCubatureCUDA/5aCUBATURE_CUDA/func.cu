
DF func Fidx[2];

DF double f0(const double *y, int index, int id)
{  
  double g;
  g=1;
          for(int i=0;i<1;i++) 
  g*=exp(-y[i]);
 
  return g;

};

DF double f1(const double *y, int index, int id)
{  
  double g;
  g=1.0;
          for(int i=0;i<2;i++) 
  g*=exp(-y[i]);
 
  return g;

};



HDF double F(const double *y, int index,int id) 
{


  Fidx[1]=&f0;
  Fidx[0]=&f1;
  
  return Fidx[index](y,index,id);
};


