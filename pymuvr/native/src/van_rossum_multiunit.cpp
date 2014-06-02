#include<cstdlib>
#include<vector>
#include<cmath>
#include<iostream>

/* Check if we're compiling with Visual Studio */
#ifdef _MSC_VER
#include <float.h>
#endif

#include "van_rossum_multiunit.hpp"

/* Visual Studio does not have isfinite; use _finite() instead */
#ifdef _MSC_VER
bool isfinite(long double x)
{
    return _finite(x);
}
#endif

using namespace std;

void d_exp_markage(double **d_matrix,vector<vector<vector<double> > > & trains, double tau,double c)
{


  /*
    If tau is zero, we have to use a different algorithm as the
    expontial kernel becomes a delta function.
  */
  if (tau==0)
    {
      d_square_zero_tau(d_matrix, trains, c);
      return;
    }

  unsigned int big_n=trains.size();
  unsigned int big_p=trains.front().size();

  // make sure d_matrix is initialised to zero
  for(unsigned int n=0;n<big_n;n++){
    for(unsigned int m=0;m<big_n;m++){
      d_matrix[n][m] = 0;
    }
  }

  vector<vector<vector<long double> > > exp_ps(big_n,vector<vector<long double> >(big_p));
  vector<vector<vector<long double> > > exp_ns(big_n,vector<vector<long double> >(big_p));

  for(unsigned int n=0;n<big_n ;++n)
    for(unsigned int p=0;p<big_p ;++p)
      try{
	expage(exp_ps[n][p], exp_ns[n][p],trains[n][p],tau);  
      }
      catch (char const * e){
	throw;
      }


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

	if (d>0){
	  d_matrix[n][m] = d_matrix[m][n] = sqrt(d);
	} // else, leave it at 0.
      }
}


void d_exp_markage_rect(double **d_matrix,
			vector<vector<vector<double> > > & trains1,
			vector<vector<vector<double> > > & trains2,
			double tau,
			double c)
{

  /*
    If tau is zero, we have to use a different algorithm as the
    expontial kernel becomes a delta function.
  */
  if (tau==0)
    {
      d_rect_zero_tau(d_matrix, trains1, trains2, c);
      return;
    }

  unsigned int big_n=trains1.size();
  unsigned int big_m=trains2.size();
  unsigned int big_p=trains1.front().size();

  if (trains2.front().size() != big_p){
    throw "d_exp_markage_rect: the observations in both lists must have the same number of cells.";
  }


  vector<vector<vector<long double> > > exp_ps1(big_n,vector<vector<long double> >(big_p));
  vector<vector<vector<long double> > > exp_ns1(big_n,vector<vector<long double> >(big_p));
  vector<vector<vector<long double> > > exp_ps2(big_m,vector<vector<long double> >(big_p));
  vector<vector<vector<long double> > > exp_ns2(big_m,vector<vector<long double> >(big_p));

  for(unsigned int n=0;n<big_n ;++n)
    for(unsigned int p=0;p<big_p ;++p)
      try{
	expage(exp_ps1[n][p], exp_ns1[n][p],trains1[n][p],tau);  
      }
      catch (char const * e){
	throw;
      }
  for(unsigned int n=0;n<big_m ;++n)
    for(unsigned int p=0;p<big_p ;++p)
      try{
	expage(exp_ps2[n][p], exp_ns2[n][p],trains2[n][p],tau);  
      }
      catch (char const * e){
	throw;
      }

  vector<vector<vector<double> > > fs1(big_n,vector<vector<double> >(big_p));
  vector<vector<vector<double> > > fs2(big_m,vector<vector<double> >(big_p));
  
  for(unsigned int n=0;n<big_n ;++n)
    for(unsigned int p=0;p<big_p ;++p)
      markage(fs1[n][p], exp_ps1[n][p], exp_ns1[n][p], trains1[n][p],tau);
  for(unsigned int n=0;n<big_m ;++n)
    for(unsigned int p=0;p<big_p ;++p)
      markage(fs2[n][p], exp_ps2[n][p], exp_ns2[n][p], trains2[n][p],tau);

  for(unsigned int n=0;n<big_n ;++n)
    for(unsigned int m=0;m<big_m ;++m)
      {
	d_matrix[n][m]=0;
	double d = 0;
	
	//R_p
	for(unsigned int p=0;p<big_p ;++p)
	  {
	    d+=big_r_with_exp_markage(fs1[n][p]);
	    d+=big_r_with_exp_markage(fs2[m][p]);
	    d-=2*big_r_with_exp_markage(trains1[n][p],fs1[n][p],exp_ps1[n][p],exp_ns1[n][p],trains2[m][p],fs2[m][p],exp_ps2[m][p],exp_ns2[m][p]);
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
	
	d+=2*c*r_pq;

	if (d>0){
	  d_matrix[n][m] = sqrt(d);
	} // else, leave it at 0.
      }
}


void d_square_zero_tau(double **d_matrix,
		       vector<vector<vector<double> > > & trains,
		       double c)
{

  unsigned int big_n=trains.size();

  for(unsigned int n=0; n<big_n; ++n)
    {
      /* Set the diagonal terms of the distance matrix to zero */
      d_matrix[n][n] = 0;
      
      /* Compute off-diagonal terms */
      for(unsigned int m=n+1; m<big_n; ++m)
	{
	  d_matrix[n][m] = d_single_observation_zero_tau(trains[n],
							 trains[m],
							 c);
	  /* The matrix must be symmetric */
	  d_matrix[m][n] = d_matrix[n][m];
	}
    }
}

void d_rect_zero_tau(double **d_matrix,
		     vector<vector<vector<double> > > & trains1,
		     vector<vector<vector<double> > > & trains2,
		     double c)
{
  unsigned int big_n=trains1.size();
  unsigned int big_m=trains2.size();

  for(unsigned int n=0; n<big_n; ++n)
    for(unsigned int m=0; m<big_m; ++m)
      {
	d_matrix[n][m] = d_single_observation_zero_tau(trains1[n],
						       trains2[m],
						       c);
      }
}



double d_single_observation_zero_tau(vector<vector<double> > & obs1,
				     vector<vector<double> > & obs2,
				     double c)
{

  unsigned int cells = obs1.size();
  if (obs2.size()!=cells)
    throw "Trying to compare two observations with a different number of cells.";

  double d_same = 0;
  double d_cross = 0;
  
  for(unsigned int p=0; p<cells; ++p)
    {
      /* Same-unit ("labelled-line") terms */
      d_same += single_unit_square_norm_zero_tau(obs1[p]);
      d_same += single_unit_square_norm_zero_tau(obs2[p]);
      d_same -= 2 * single_unit_scalar_product_zero_tau(obs1[p],
							obs2[p]);
      
      /* Cross-unit ("summed-population") terms */
      for(unsigned int q=p+1; q<cells; ++q)
	{
	  /* Same-observation, cross-unit */
	  d_cross += single_unit_scalar_product_zero_tau(obs1[p],
							 obs1[q]);
	  d_cross += single_unit_scalar_product_zero_tau(obs2[p],
							 obs2[q]);
	  /* Cross-observation, cross-unit */
	  d_cross -= single_unit_scalar_product_zero_tau(obs1[p],
							 obs2[q]);
	  d_cross -= single_unit_scalar_product_zero_tau(obs1[q],
							 obs2[p]);
	}
    }
  
  /*
    Sum same-unit and cross-unit with the appropriate weighting
    to give the desired interpolation between labelled-line and
    summed-population distance.
  */
  double d = d_same + 2 * c * d_cross;

  if (d>0)
    {
      return sqrt(d);
    } else
    {
      return 0;
    }
}



double big_r_with_exp_markage(vector<double> & fs)
{
  unsigned int f_size=fs.size();

  if(f_size==0)
    return 0;

  double norm=f_size;
  
  for(unsigned int i=1;i<f_size ;++i)
    norm+=2*fs[i];
   
  return norm;

}


int single_unit_square_norm_zero_tau(vector<double> & train)
{
  return train.size();
}


int single_unit_scalar_product_zero_tau(vector<double> & train_a,
					vector<double> & train_b)
{

  unsigned int train_a_size=train_a.size();
  unsigned int train_b_size=train_b.size();

  if(train_a_size==0 || train_b_size==0)
    return 0;

  int count = 0;
  int place = train_a_size - 1;

  for (int i=train_b_size-1; i>=0; --i)
    {
      /*
	Look for the index of largest spike time in train_a which is
	smaller or equal than the ith spike time of train_b. Leave the
	index set to zero if you don't find any.
      */
      while(place>0 && train_a[place]>train_b[i])
	place--;
      /*
	If the spike selected in train_a coincides with the one we're
	considering in train_b, add 1 to the pair count.
      */
      if (train_a[place]==train_b[i])
	{
	  count+=1;
	}
    }

  return count;
}


double big_r_with_exp_markage(vector<double> & train_a, vector<double> & f_a,vector<long double> & e_pos_a,vector<long double> & e_neg_a,vector<double> & train_b, vector<double> & f_b,vector<long double> & e_pos_b,vector<long double> & e_neg_b)
{
  
  unsigned int train_a_size=train_a.size();
  unsigned int train_b_size=train_b.size();

  if(train_a_size==0||train_b_size==0)
    return 0;


  double x=0;

  int place=train_a_size-1;

  for(int i=train_b_size-1;i>=0 ;--i)
    {
      /*
	Look for the index of largest spike time in train_a which is
	smaller or equal than the ith spike time of train_b. Exit from
	the loop if you don't find any. Note that this, in general, is
	not the same thing as calculating J(i) in the paper (Houghton
	and Kreuz 2012), because we allow for train_a[place] to be
	equal to train_b[i].
      */
      while(place>=0&&train_a[place]>train_b[i])
	  place--;
      if(place<0)
	break;
      /*
	If the spike selected in train_a coincides with the one we're
	considering in train_b, add 1 to the quantity being
	calculated. Adjust the 'place' index to point at the previous
	spike, which now is guaranteed to correspond to the J(i) index
	in the paper mentioned above; unless it does not exist, in
	which case we exit from the loop. Note that this passage is
	not mirrored in the 'symmetric' loop below, as we only need to
	count each pair of coincident spikes once.
      */
      if (train_a[place]==train_b[i])
	{
	  x+=1;
	  place--;
	}
      if(place<0)
	break;
      /*
	Compute the ith term of the main sum we're calculating, using
	the precomputed exponentials and the markage vector.
      */
      x+=e_pos_a[place]*e_neg_b[i]*(f_a[place]+1);	
    }

  place=train_b_size-1;

  for(int i=train_a_size-1;i>=0 ;--i)
    {
      /*
	Unlike above, calculate directly the index corresponding to
	the largest spike time in train_b which is _strictly_ smaller
	than the ith spike time of train_a. Exit from the loop if you
	don't find any.
      */
      while(place>=0&&train_b[place]>=train_a[i])
	place--;
      if(place<0)
	break;
      /*
	Compute the ith term of the main sum we're calculating, using
	the precomputed exponentials and the markage vector.
      */
      x+=e_pos_b[place]*e_neg_a[i]*(f_b[place]+1);
    }
  return x;

}



void expage(vector<long double> & e_pos, vector<long double> & e_neg,vector<double> & train,double tau)
{

  unsigned int train_size=train.size();

  e_pos.resize(train_size);
  e_neg.resize(train_size);
  
  for(unsigned int i=0;i<train_size ;++i)
    {
      e_pos[i]=exp((long double)(train[i]/tau));
      e_neg[i]=exp((long double)(-train[i]/tau));
    }
  
  /*Check if any over/underflow occurred while calculating the
    exponentials*/
  if (train_size>0 && (e_neg.front()==0 || !isfinite(e_pos.back())))
    {
      throw "tau is too small compared to the spike times. Please use a larger value for tau, or shorter spike trains.";
    }

}

void markage(vector<double> & f, vector<long double> & e_pos, vector<long double> & e_neg, vector<double> & train,double tau)
{

  unsigned int train_size=train.size();
  
  f.resize(train_size);
  
  if (train_size == 0)
    return;

  f[0]=0;

  for(unsigned int i=1;i<train_size ;++i)
    f[i]=(1+f[i-1])*e_pos[i-1]*e_neg[i];

}
