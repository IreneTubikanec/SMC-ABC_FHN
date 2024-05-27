#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;
using namespace Eigen;

//#Matrix-Vector Multiplication
// [[Rcpp::export]]
NumericVector mv_multFHN_(NumericMatrix mat, NumericVector vec)
{
  NumericVector ret(mat.nrow());
  double temp=0;
  
  for(int i = 0; i < mat.nrow(); i++)
  {
    for(int j = 0; j < vec.size(); j++)
    {
      temp = temp + mat(i,j) * vec[j];
    }
    ret[i]=temp;
    temp=0;
  }
  return ret;
};

//linear function (of splitting method)
// [[Rcpp::export]]
NumericVector linearFHN_Cpp_(NumericVector vec, NumericMatrix dm, NumericMatrix cm, NumericVector randvec)
{
  NumericVector ret=mv_multFHN_(dm,vec)+mv_multFHN_(cm,randvec);
  return ret;
};

//nonlinear function (of splitting method)
// [[Rcpp::export]]
NumericVector nonlinearFHN_Cpp_(NumericVector vec, double h, double eps, double beta)
{
  double a=1.0/eps;
  NumericVector ret(2);
  ret(0)=vec(0)/sqrt((exp(-2.0*a*h)+vec(0)*vec(0)*(1.0-exp(-2.0*a*h))));
  ret(1)=vec(1)+beta*h;
  return ret;
};

//#Strang splitting method
// [[Rcpp::export]]
NumericMatrix FHN_Splitting_Cpp_(NumericVector grid_i, double h_i, NumericVector startv_i, NumericMatrix dm_i, NumericMatrix cm_i,  double eps_i, double beta_i)
{
  double h=h_i;
  NumericVector startv=startv_i;
  NumericVector grid=grid_i;
  int iter=grid.size();
  
  NumericMatrix randarr(2,iter);
  randarr(0,_)=rnorm(iter);
  randarr(1,_)=rnorm(iter);
  
  NumericMatrix dm=dm_i;
  NumericMatrix cm=cm_i;
  
  double eps=eps_i;
  double beta=beta_i;
  
  NumericMatrix sol(2,iter);
  sol(_, 0)=startv;
  NumericVector newv=startv;
  NumericVector randvec(2);
  
  for(int i=1;i<iter;i++)
  {
    randvec=randarr(_,i);
    newv=nonlinearFHN_Cpp_(newv,h/2,eps,beta);
    newv=linearFHN_Cpp_(newv,dm,cm,randvec);
    newv=nonlinearFHN_Cpp_(newv,h/2,eps,beta);
    sol(_,i)=newv;
  }
  
  NumericMatrix ret=sol;
  
  return sol;
};




