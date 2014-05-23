#ifndef VAN_ROSSUM_H
#define VAN_ROSSUM_H

#include<vector>

using namespace std;


double d(vector<vector<double> > & trains_a,vector<vector<double> > & trains_b, double tau,double c);

//trains[n][p][i]
void d(double **d_matrix,vector<vector<vector<double> > > & trains, double tau,double c);
void d_exp(double **d_matrix,vector<vector<vector<double> > > & trains, double tau,double c);
void d_exp_markage(double **d_matrix,vector<vector<vector<double> > > & trains, double tau,double c);
void d_exp_markage_rect(double **d_matrix,vector<vector<vector<double> > > & trains1,vector<vector<vector<double> > > & trains2, double tau,double c);

double big_r(vector<double> & u,double tau);
double big_r_with_exp(vector<double> & u,vector<long double> & exp_p,vector<long double> & exp_m,double tau);
double big_r_with_exp_markage(vector<double> & fs);

double big_r(vector<double> & u_a,vector<double> & u_b,double tau);
double big_r_with_exp(vector<double> & u_a,vector<long double> & exp_p_a,vector<long double> & exp_m_a,vector<double> & u_b,vector<long double> & exp_p_b,vector<long double> & exp_m_b,double tau);
double big_r_with_exp_markage(vector<double> & train_a, vector<double> & f_a,vector<long double> & e_pos_a,vector<long double> & e_neg_a,vector<double> & train_b, vector<double> & f_b,vector<long double> & e_pos_b,vector<long double> & e_neg_b);

void expage(vector<long double> & e_pos, vector<long double> & e_neg,vector<double> & train,double tau);
void markage(vector<double> & f, vector<long double> & e_pos, vector<long double> & e_neg, vector<double> & train,double tau);

#endif
