#include<cstdlib>
#include<vector>
#include<cmath>
#include<iostream>

#include "Van_Rossum_Multiunit.hpp"

using namespace std;

//trains[p][i]

//single train
double d(vector<vector<double> > & trains_a,vector<vector<double> > & trains_b, double tau,double c)
{

  double big_p=trains_a.size();

  double d=0; 

  //R_p
  for(unsigned int p=0;p<big_p ;++p)
    {
      d+=big_r(trains_a[p],tau);
      d+=big_r(trains_b[p],tau);
      d-=2*big_r(trains_a[p],trains_b[p],tau);
    }

  double r_pq=0;

  for(unsigned int p=0;p<big_p-1 ;++p)
    for(unsigned int q=p+1;q<big_p ;++q)
      {
	r_pq+=big_r(trains_a[p],trains_a[q],tau);
	r_pq+=big_r(trains_b[p],trains_b[q],tau);
	r_pq-=big_r(trains_a[p],trains_b[q],tau);
	r_pq-=big_r(trains_b[p],trains_a[q],tau);
      }

  d+=2*c*r_pq;

  return sqrt(d);

}

//trains[n][p][i]

void d(double **d_matrix,vector<vector<vector<double> > > & trains, double tau,double c)
{
  int big_n=trains.size();
  int big_p=trains.front().size();

  vector<double> squares(big_n,0.0);

  for(unsigned int n=0;n<big_n;++n)
    {
      for(unsigned int p=0;p<big_p ;++p)
	squares[n]+=big_r(trains[n][p],tau);
      for(unsigned int p=0;p<big_p-1 ;++p)
	for(unsigned int q=p+1;q<big_p ;++q)
	  squares[n]+=2*c*big_r(trains[n][p],trains[n][q],tau);
    }

  for(unsigned int n=0;n<big_n-1 ;++n)
    for(unsigned int m=n+1;m<big_n ;++m)
      {
	double d=squares[n]+squares[m];
	for(unsigned int p=0;p<big_p ;++p)
	  d-=2*big_r(trains[n][p],trains[m][p],tau);
	  
	for(unsigned int p=0;p<big_p-1 ;++p)
	  for(unsigned int q=p+1;q<big_p ;++q)
	    {
	      d-=2*c*big_r(trains[n][p],trains[m][q],tau);
	      d-=2*c*big_r(trains[m][p],trains[n][q],tau);
	    }

	d_matrix[n][m]=sqrt(2*d/tau);
	d_matrix[m][n]=sqrt(2*d/tau);

      }
}


//trains[n][p][i]

void d_exp(double **d_matrix,vector<vector<vector<double> > > & trains, double tau,double c)
{

  int big_n=trains.size();
  int big_p=trains.front().size();

  vector<vector<vector<double> > > exp_ps(big_n,vector<vector<double> >(big_p));
  vector<vector<vector<double> > > exp_ns(big_n,vector<vector<double> >(big_p));
 
  for(unsigned int n=0;n<big_n ;++n)
    for(unsigned int p=0;p<big_p ;++p)
      expage(exp_ps[n][p], exp_ns[n][p],trains[n][p],tau);  

  vector<double> squares(big_n,0.0);

  for(unsigned int n=0;n<big_n;++n)
    {
      for(unsigned int p=0;p<big_p ;++p)
	squares[n]+=big_r_with_exp(trains[n][p],exp_ps[n][p],exp_ns[n][p],tau);
      for(unsigned int p=0;p<big_p-1 ;++p)
	for(unsigned int q=p+1;q<big_p ;++q)
	  squares[n]+=2*c*big_r_with_exp(trains[n][p],exp_ps[n][p],exp_ns[n][p],trains[n][q],exp_ps[n][q],exp_ns[n][q],tau);
    }

  for(unsigned int n=0;n<big_n-1 ;++n)
    for(unsigned int m=n+1;m<big_n ;++m)
      {
	double d=squares[n]+squares[m];
	for(unsigned int p=0;p<big_p ;++p)
	  d-=2*big_r_with_exp(trains[n][p],exp_ps[n][p],exp_ns[n][p],trains[m][p],exp_ps[m][p],exp_ns[m][p],tau);
	  
	for(unsigned int p=0;p<big_p-1 ;++p)
	  for(unsigned int q=p+1;q<big_p ;++q)
	    {
	      d-=2*c*big_r_with_exp(trains[n][p],exp_ps[n][p],exp_ns[n][p],trains[m][q],exp_ps[m][q],exp_ns[m][q],tau);
	      d-=2*c*big_r_with_exp(trains[m][p],exp_ps[m][p],exp_ns[m][p],trains[n][q],exp_ps[n][q],exp_ns[n][q],tau);
	    }

	d_matrix[n][m]=sqrt(d);
	d_matrix[m][n]=sqrt(d);

      }
}


//trains[n][p][i]

void d_exp_markage(double **d_matrix,vector<vector<vector<double> > > & trains, double tau,double c)
{

  int big_n=trains.size();
  int big_p=trains.front().size();

  // make sure d_matrix is initialised to zero
  for(unsigned int n=0;n<big_n;n++){
    for(unsigned int m=0;m<big_n;m++){
      d_matrix[n][m] = 0;
    }
  }

  vector<vector<vector<double> > > exp_ps(big_n,vector<vector<double> >(big_p));
  vector<vector<vector<double> > > exp_ns(big_n,vector<vector<double> >(big_p));

  for(unsigned int n=0;n<big_n ;++n)
    for(unsigned int p=0;p<big_p ;++p)
      expage(exp_ps[n][p], exp_ns[n][p],trains[n][p],tau);  

  vector<vector<vector<double> > > fs(big_n,vector<vector<double> >(big_p));
  
  for(unsigned int n=0;n<big_n ;++n)
    for(unsigned int p=0;p<big_p ;++p)
      markage(fs[n][p], exp_ps[n][p], exp_ns[n][p], trains[n][p],tau);

  vector<double> squares(big_n,0.0);

  for(unsigned int n=0;n<big_n;++n)
    {
      for(unsigned int p=0;p<big_p ;++p)
	squares[n]+=big_r_with_exp_markage(fs[n][p]);
      for(unsigned int p=0;p<big_p-1 ;++p)
	for(unsigned int q=p+1;q<big_p ;++q)
	  squares[n]+=2*c*big_r_with_exp_markage(trains[n][p],fs[n][p],exp_ps[n][p],exp_ns[n][p],trains[n][q],fs[n][q],exp_ps[n][q],exp_ns[n][q]);
    }

  for(unsigned int n=0;n<big_n-1 ;++n)
    for(unsigned int m=n+1;m<big_n ;++m)
      {
	double d=squares[n]+squares[m];
	for(unsigned int p=0;p<big_p ;++p)
	  d-=2*big_r_with_exp_markage(trains[n][p],fs[n][p],exp_ps[n][p],exp_ns[n][p],trains[m][p],fs[m][p],exp_ps[m][p],exp_ns[m][p]);
	 	  
	for(unsigned int p=0;p<big_p-1 ;++p)
	  for(unsigned int q=p+1;q<big_p ;++q)
	    {
	      d-=2*c*big_r_with_exp_markage(trains[n][p],fs[n][p],exp_ps[n][p],exp_ns[n][p],trains[m][q],fs[m][q],exp_ps[m][q],exp_ns[m][q]);
	      d-=2*c*big_r_with_exp_markage(trains[m][p],fs[m][p],exp_ps[m][p],exp_ns[m][p],trains[n][q],fs[n][q],exp_ps[n][q],exp_ns[n][q]);
	    }

	d_matrix[n][m]=sqrt(d);
	d_matrix[m][n]=sqrt(d);

      }
}

void d_exp_markage_rect(double **d_matrix,
			vector<vector<vector<double> > > & trains1,
			vector<vector<vector<double> > > & trains2,
			double tau,
			double c)
{

  int big_n=trains1.size();
  int big_m=trains2.size();
  int big_p=trains1.front().size();

  if (trains2.front().size() != big_p){
    cout << "d_exp_markage_rect: the observations in both lists must have the same number of cells." << endl;
  }


  vector<vector<vector<double> > > exp_ps1(big_n,vector<vector<double> >(big_p));
  vector<vector<vector<double> > > exp_ns1(big_n,vector<vector<double> >(big_p));
  vector<vector<vector<double> > > exp_ps2(big_m,vector<vector<double> >(big_p));
  vector<vector<vector<double> > > exp_ns2(big_m,vector<vector<double> >(big_p));

  for(unsigned int n=0;n<big_n ;++n)
    for(unsigned int p=0;p<big_p ;++p)
      expage(exp_ps1[n][p], exp_ns1[n][p],trains1[n][p],tau);  
  for(unsigned int n=0;n<big_m ;++n)
    for(unsigned int p=0;p<big_p ;++p)
      expage(exp_ps2[n][p], exp_ns2[n][p],trains2[n][p],tau);  

  vector<vector<vector<double> > > fs1(big_n,vector<vector<double> >(big_p));
  vector<vector<vector<double> > > fs2(big_m,vector<vector<double> >(big_p));
  
  for(unsigned int n=0;n<big_n ;++n)
    for(unsigned int p=0;p<big_p ;++p)
      markage(fs1[n][p], exp_ps1[n][p], exp_ns1[n][p], trains1[n][p],tau);
  for(unsigned int n=0;n<big_m ;++n)
    for(unsigned int p=0;p<big_p ;++p)
      markage(fs2[n][p], exp_ps2[n][p], exp_ns2[n][p], trains2[n][p],tau);

  for(unsigned int n=0;n<big_n-1 ;++n)
    for(unsigned int m=0;m<big_m-1 ;++m)
      {
	d_matrix[n][m]=0; 
	
	//R_p
	for(unsigned int p=0;p<big_p ;++p)
	  {
	    d_matrix[n][m]+=big_r_with_exp_markage(fs1[n][p]);
	    d_matrix[n][m]+=big_r_with_exp_markage(fs2[m][p]);
	    d_matrix[n][m]-=2*big_r_with_exp_markage(trains1[n][p],fs1[n][p],exp_ps1[n][p],exp_ns1[n][p],trains2[m][p],fs2[m][p],exp_ps2[m][p],exp_ns2[m][p]);
	  }
	
	double r_pq=0;
	
	for(unsigned int p=0;p<big_p-1 ;++p)
	  for(unsigned int q=p+1;q<big_p ;++q)
	    {
	      r_pq+=big_r_with_exp_markage(trains1[n][p],fs1[n][p],exp_ps1[n][p],exp_ns1[n][p],trains1[n][q],fs1[n][q],exp_ps1[n][q],exp_ns1[n][q]);
	      r_pq+=big_r_with_exp_markage(trains2[m][p],fs2[m][p],exp_ps2[m][p],exp_ns2[m][p],trains2[m][q],fs2[m][q],exp_ps2[m][q],exp_ns2[m][q]);
	      r_pq-=big_r_with_exp_markage(trains1[n][p],fs1[n][p],exp_ps1[n][p],exp_ns1[n][p],trains2[m][q],fs2[m][q],exp_ps2[m][q],exp_ns2[m][q]);
	      r_pq-=big_r_with_exp_markage(trains2[m][p],fs2[m][p],exp_ps2[m][p],exp_ns2[m][p],trains1[n][q],fs1[n][q],exp_ps1[n][q],exp_ns1[n][q]);
	    }
	
	d_matrix[n][m]+=2*c*r_pq;
	d_matrix[n][m]=sqrt(d_matrix[n][m]);
      }
}


double big_r(vector<double> & u,double tau)
{
  double u_size=u.size();

  if(u_size==0)
    return 0;

  double total=u_size;

  for(unsigned int i=0;i<u_size ;++i)
    for(unsigned int j=i+1;j<u_size; j++)
      total+=2*exp(-(u[j]-u[i])/tau);

  return total;

}


double big_r_with_exp(vector<double> & u,vector<double> & exp_p,vector<double> & exp_m,double tau)
{
  double u_size=u.size();

  if(u_size==0)
    return 0;

  double total=u_size;

  for(unsigned int i=0;i<u_size ;++i)
    for(unsigned int j=i+1;j<u_size; j++)
      total+=2*exp_p[i]*exp_m[j];

  return total;

}


double big_r_with_exp_markage(vector<double> & fs)
{
  int f_size=fs.size();

  if(f_size==0)
    return 0;

  double norm=f_size;
  
  for(unsigned int i=1;i<f_size ;++i)
    norm+=2*fs[i];
   
  return norm;

}


double big_r(vector<double> & u_a,vector<double> & u_b,double tau)
{
  int u_a_size=u_a.size();
  int u_b_size=u_b.size();


  if(u_a_size==0||u_b_size==0)
    return 0;

  double total=0;

  for(unsigned int i=0;i<u_a_size ;++i)
    for(unsigned int j=0;j<u_b_size; j++)
      total+=exp(-fabs(u_a[i]-u_b[j])/tau);

  return total;

}


double big_r_with_exp(vector<double> & u_a,vector<double> & exp_p_a,vector<double> & exp_m_a,vector<double> & u_b,vector<double> & exp_p_b,vector<double> & exp_m_b,double tau)
{
  int u_a_size=u_a.size();
  int u_b_size=u_b.size();

  if(u_a_size==0||u_b_size==0)
    return 0;

  double total=0;

  for(unsigned int i=0;i<u_a_size ;++i)
    for(unsigned int j=0;j<u_b_size; j++)
      if(u_a[i]>u_b[j])
	total+=exp_m_a[i]*exp_p_b[j];
      else
	total+=exp_p_a[i]*exp_m_b[j];

  return total;

}

double big_r_with_exp_markage(vector<double> & train_a, vector<double> & f_a,vector<double> & e_pos_a,vector<double> & e_neg_a,vector<double> & train_b, vector<double> & f_b,vector<double> & e_pos_b,vector<double> & e_neg_b)
{
  
  int train_a_size=train_a.size();
  int train_b_size=train_b.size();

  if(train_a_size==0||train_b_size==0)
    return 0;


  double x=0;

  int place=train_a_size-1;

  for(int i=train_b_size-1;i>=0 ;--i)
    {
      while(place>=0&&train_a[place]>train_b[i])
	  place--;
      if(place<0)
	break;
      x+=e_pos_a[place]*e_neg_b[i]*(f_a[place]+1);	
    }

  place=train_b_size-1;

  for(int i=train_a_size-1;i>=0 ;--i)
    {
      while(place>=0&&train_b[place]>train_a[i])
	place--;
      if(place<0)
	break;
      x+=e_pos_b[place]*e_neg_a[i]*(f_b[place]+1);
    }

  return x;

}



void expage(vector<double> & e_pos, vector<double> & e_neg,vector<double> & train,double tau)
{

  int train_size=train.size();

  e_pos.resize(train_size);
  e_neg.resize(train_size);
  
  for(unsigned int i=0;i<train_size ;++i)
    {
      e_pos[i]=(exp(train[i]/tau));
      e_neg[i]=(exp(-train[i]/tau));
    }
  

}

void markage(vector<double> & f, vector<double> & e_pos, vector<double> & e_neg, vector<double> & train,double tau)
{

  int train_size=train.size();
  
  f.resize(train_size);
  
  f[0]=0;

  for(unsigned int i=1;i<train_size ;++i)
    f[i]=(1+f[i-1])*e_pos[i-1]*e_neg[i];

}
