#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <map>
#include "GAIT.h"
#include "MyMath.h"
#include "CenData.h"
#include "Region.h"
#include "ExpData.h"
#include "OPT.h"
#include "PHI.h"
#include "student_t_distribution.h"


using namespace std;


GAIT::GAIT()
{
	N_tgen=0;
}
GAIT::~GAIT()
{

}

double opt2dArea(map<Region,deninfo_t,RegionCompare> *S,double X1,double X2, double Y1, double Y2,double max1, double max2)
{
	double a1=0;
	double a2=0;
	double b1=0;
	double b2=0;
	double x;
	double y;
	double P=0;
	double L=0;
	double Den=0;
	if(max1>=max2)	L=10*max1;
	if(max1<max2)	L=10*max2;


	if(X2!=INFINITY&&Y2!=INFINITY)
	{
		for( map<Region,deninfo_t,RegionCompare>::iterator it=S->begin() ; it!=S->end() ; it++ )
		{
			a1=it->first.at(0).begin;
			a2=it->first.at(0).end;
			b1=it->first.at(1).begin;
			b2=it->first.at(1).end;
			if(a1<=X1)	a1=X1;
			if(a2>=X2)	a2=X2;
			if(b1<=Y1)	b1=Y1;
			if(b2>=Y2)	b2=Y2;
			x=a2-a1;
			y=b2-b1;
			if(x<0)	x=0;
			if(y<0)	y=0;
			P=P+x*y*it->second.den;
		}
		return P;
	}
	if(X2==INFINITY&&Y2!=INFINITY)
	{
		X2=L;
		for( map<Region,deninfo_t,RegionCompare>::iterator it=S->begin() ; it!=S->end() ; it++ )
		{
			a1=it->first.at(0).begin;
			a2=it->first.at(0).end;
			b1=it->first.at(1).begin;
			b2=it->first.at(1).end;

			if(it->first.at(0).end==INFINITY)
			{
				a2=L;
			}
			Den=it->second.prob/(a2-a1)/(b2-b1);
			if(a1<=X1)	a1=X1;
			if(a2>=X2)	a2=X2;
			if(b1<=Y1)	b1=Y1;
			if(b2>=Y2)	b2=Y2;
			x=a2-a1;
			y=b2-b1;
			if(x<0)	x=0;
			if(y<0)	y=0;
			P=P+x*y*Den;
		}
		return P;

	}
	if(X2!=INFINITY&&Y2==INFINITY)
	{
		Y2=L;
		for( map<Region,deninfo_t,RegionCompare>::iterator it=S->begin() ; it!=S->end() ; it++ )
		{
			a1=it->first.at(0).begin;
			a2=it->first.at(0).end;
			b1=it->first.at(1).begin;
			b2=it->first.at(1).end;

			if(it->first.at(1).end==INFINITY)
			{
				b2=L;
			}
			Den=it->second.prob/(a2-a1)/(b2-b1);
			if(a1<=X1)	a1=X1;
			if(a2>=X2)	a2=X2;
			if(b1<=Y1)	b1=Y1;
			if(b2>=Y2)	b2=Y2;
			x=a2-a1;
			y=b2-b1;
			if(x<0)	x=0;
			if(y<0)	y=0;
			P=P+x*y*Den;
		}
		return P;

	}
	X2=L;
	Y2=L;
	for( map<Region,deninfo_t,RegionCompare>::iterator it=S->begin() ; it!=S->end() ; it++ )
	{
		a1=it->first.at(0).begin;
		a2=it->first.at(0).end;
		b1=it->first.at(1).begin;
		b2=it->first.at(1).end;
		if(it->first.at(0).end==INFINITY)
		{
			a2=L;
		}

		if(it->first.at(1).end==INFINITY)
		{
			b2=L;
		}
		Den=it->second.prob/(a2-a1)/(b2-b1);
		if(a1<=X1)	a1=X1;
		if(a2>=X2)	a2=X2;
		if(b1<=Y1)	b1=Y1;
		if(b2>=Y2)	b2=Y2;
		x=a2-a1;
		y=b2-b1;
		if(x<0)	x=0;
		if(y<0)	y=0;
		P=P+x*y*Den;
	}

	return P;
}


double unifRand()
{
    return rand() / double(RAND_MAX);
}

double unifRand(double a, double b)
{
    return (b-a)*unifRand() + a;
}

long unifRand(long n)
{

    if (n < 0) n = -n;
    if (n==0) return 0;
    long guard = (long) (unifRand() * n) +1;
    return (guard > n)? n : guard;
}

void seed()
{
    srand(time(NULL));
}

double rand_val(int seed)
{
  const long  a =      16807;  // Multiplier
  const long  m = 2147483647;  // Modulus
  const long  q =     127773;  // m div a
  const long  r =       2836;  // m mod a
  static long x;               // Random int value
  long        x_div_q;         // x divided by q
  long        x_mod_q;         // x modulo q
  long        x_new;           // New x value

  // Set the seed if argument is non-zero and then return zero
  if (seed > 0)
  {
    x = seed;
    return(0.0);
  }

  // RNG using integer arithmetic
  x_div_q = x / q;
  x_mod_q = x % q;
  x_new = (a * x_mod_q) - (r * x_div_q);
  if (x_new > 0)
    x = x_new;
  else
    x = x_new + m;

  // Return a random value between 0.0 and 1.0
  return((double) x / m);
}

double expon(double x)
{
  double z;                     // Uniform random number (0 < z < 1)
  double exp_value;             // Computed exponential value to be returned

  // Pull a uniform random number (0 < z < 1)
  do
  {
    z = rand_val(0);
  }
  while ((z == 0) || (z == 1));

  // Compute exponential random variable using inversion method
  exp_value = -x * log(z);

  return(exp_value);
}



void GAIT::Tgen(int NDIV,int N_sample,CenData* Data,Density* Den, bool Bound)
{


	double S_min_x=INFINITY;
	double S_min_y=INFINITY;
	double S_max_x=-INFINITY;
	double S_max_y=-INFINITY;
	double dim1;
	double dim2;

	vector<double> bin(NDIV+1,0);
	vector<double> bin2(NDIV+1,0);
	vector<double> p(NDIV,0);
	vector<double> p2(NDIV,0);
	vector<double> cp(NDIV+1,0);
	vector<double> cp2(NDIV+1,0);
	vector<double> x(2*N_sample,0);
	vector<double> t(2*N_sample,0);
	double Psum=0;
	vector<double> Y(NDIV,0);
	vector<vector<double> >	pp(NDIV,Y);

	vector<cendata_t> *inputdata_T1 = Data->getData(0);
	vector<cendata_t> *inputdata_T2 = Data->getData(1);

	map<Region,deninfo_t,RegionCompare> *Sdata = Den->getDensity();
	for( map<Region,deninfo_t,RegionCompare>::iterator it=Sdata->begin() ; it!=Sdata->end() ; it++ )
	{
		if(it->first.at(0).begin<S_min_x){S_min_x=it->first.at(0).begin;}
		if(it->first.at(1).begin<S_min_y){S_min_y=it->first.at(1).begin;}
		if(it->first.at(0).begin>S_max_x){S_max_x=it->first.at(0).begin;}
		if(it->first.at(1).begin>S_max_y){S_max_y=it->first.at(1).begin;}
	}
	dim1=(S_max_x-S_min_x)/2/(NDIV-1);
	dim2=(S_max_y-S_min_y)/2/(NDIV-1);
	vector<double> v(inputdata_T1->size(),0);
	vector<vector<double> >	TGen1(N_sample,v);
	vector<vector<double> >	TGen2(N_sample,v);

	for( unsigned int j=0 ; j<inputdata_T1->size() ; j++ )
	{
		tgen_t tgen;
		seed();
		Psum=0;
		for(int i=0;i<NDIV;i++)
		{
			p[i]=0;
			p2[i]=0;
			cp[i]=0;
			cp2[i]=0;
		}
		cp[NDIV]=0;
		cp2[NDIV]=0;
		if(inputdata_T1->at(j).e==1 && inputdata_T2->at(j).e==1)
		{
			for(int i=0; i<N_sample; i++)
			{
				tgen.T1Gen.push_back(inputdata_T1->at(j).x);
				tgen.T2Gen.push_back(inputdata_T2->at(j).x);
			}
		}
		else if(inputdata_T1->at(j).e==0 && inputdata_T2->at(j).e==1)
		{
			for(unsigned int k=0; k<bin.size()-1; k++)
			{
				bin[k]=inputdata_T1->at(j).x + (S_max_x - inputdata_T1->at(j).x)/(NDIV-1)*k;
			}
			bin[NDIV]=INFINITY;

			for(int i=0;i<NDIV;i++)
			{
				p[i]=opt2dArea(Sdata,bin[i],bin[i+1],inputdata_T2->at(j).x-dim2,inputdata_T2->at(j).x+dim2,S_max_x,S_max_y);
				Psum+=p[i];
			}



			for(int i=0;i<NDIV;i++)
			{
				p[i]=p[i]/Psum;
				for(int k=0;k<=i;k++)
				{
					cp[i+1]+=p[k];
				}
			}




			for(int i=0;i<N_sample;i++)
			{
				x[i]=t[i]=unifRand();
			}

			for(int i=0;i<NDIV;i++)
			{
				for(int k=0;k<N_sample;k++)
				{
					if(cp[i]<=x[k])
					{
						if(x[k]<cp[i+1])
						{
							if(bin[i+1]!=INFINITY)
							{
								t[k]=(bin[i+1]-bin[i])/(cp[i+1]-cp[i])*(x[k]-cp[i])+bin[i];
							}
							else
							{
								t[k]=bin[i];
							}
						}
					}
				}

			}

			for(int i=0;i<N_sample;i++)
			{
				tgen.T1Gen.push_back(t[i]);
				tgen.T2Gen.push_back(inputdata_T2->at(j).x);
			}

		}
		else if(inputdata_T1->at(j).e==1 && inputdata_T2->at(j).e==0)
		{
			for(unsigned int k=0; k<bin.size()-1; k++)
			{
				bin[k]=inputdata_T2->at(j).x + (S_max_y - inputdata_T2->at(j).x)/(NDIV-1)*k;
			}
			bin[NDIV]=INFINITY;


			for(int i=0;i<NDIV;i++)
			{
				p[i]=opt2dArea(Sdata,inputdata_T1->at(j).x-dim1,inputdata_T1->at(j).x+dim1,bin[i],bin[i+1],S_max_x,S_max_y);
				Psum+=p[i];
			}

			for(int i=0;i<NDIV;i++)
			{
				p[i]=p[i]/Psum;
				for(int k=0;k<=i;k++)
				{
					cp[i+1]+=p[k];
				}
			}
			for(int i=0;i<N_sample;i++)
			{
				x[i]=t[i]=unifRand();
			}

			for(int i=0;i<NDIV;i++)
			{
				for(int k=0;k<N_sample;k++)
				{
					if(cp[i]<=x[k])
					{
						if(x[k]<cp[i+1])
						{
							if(bin[i+1]!=INFINITY)
							{
								t[k]=(bin[i+1]-bin[i])/(cp[i+1]-cp[i])*(x[k]-cp[i])+bin[i];
							}
							else
							{
								t[k]=bin[i];
							}
						}
					}
				}

			}

			for(int i=0;i<N_sample;i++)
			{
				tgen.T1Gen.push_back(inputdata_T1->at(j).x);
				tgen.T2Gen.push_back(t[i]);
			}
		}
		else if(inputdata_T1->at(j).e==0 && inputdata_T2->at(j).e==0)
		{
			for(unsigned int k=0; k<bin.size()-1; k++)
			{
				bin[k]=inputdata_T1->at(j).x + (S_max_x - inputdata_T1->at(j).x)/(NDIV-1)*k;
				bin2[k]=inputdata_T2->at(j).x + (S_max_y - inputdata_T2->at(j).x)/(NDIV-1)*k;
			}
			bin[NDIV]=INFINITY;
			bin2[NDIV]=INFINITY;


			for(int i=0;i<NDIV;i++)
			{
				for(int ii=0;ii<NDIV;ii++)
				{
					pp[i][ii]=opt2dArea(Sdata,bin[i],bin[i+1],bin2[ii],bin2[ii+1],S_max_x,S_max_y);
					Psum+=pp[i][ii];
				}
			}

			for(int i=0;i<NDIV;i++)
			{
				for(int ii=0;ii<NDIV;ii++)
				{
					pp[i][ii]=pp[i][ii]/Psum;
				}
			}

			for(int i=0;i<NDIV;i++)
			{
				for(int ii=0;ii<NDIV;ii++)
				{
					p[i]=p[i]+pp[i][ii];
					p2[i]=p2[i]+pp[ii][i];
				}
				for(int k=0;k<=i;k++)
				{
					cp[i+1]+=p[k];
					cp2[i+1]+=p2[k];
				}
			}

			for(int i=0;i<2*N_sample;i++)
			{
				x[i]=t[i]=unifRand();
			}

			for(int i=0;i<NDIV;i++)
			{
				for(int ii=0;ii<NDIV;ii++)
				{
					for(int k=0;k<N_sample;k++)
					{
						if(cp[i]<=x[k])
						{
							if(x[k]<cp[i+1])
							{
								if(cp2[ii]<=x[N_sample+k])
								{
									if(x[N_sample+k]<cp2[ii+1])
									{
										if(bin[i+1]!=INFINITY)
										{
											t[k]=(bin[i+1]-bin[i])/(cp[i+1]-cp[i])*(x[k]-cp[i])+bin[i];
										}
										else
										{
											t[k]=bin[i];
										}
										if(bin2[ii+1]!=INFINITY)
										{
											t[N_sample+k]=(bin2[ii+1]-bin2[ii])/(cp2[ii+1]-cp2[ii])*(x[N_sample+k]-cp2[ii])+bin2[ii];
										}
										else
										{
											t[N_sample+k]=bin2[ii];
										}
									}
								}

							}
						}
					}

				}
			}

			for(int i=0;i<N_sample;i++)
			{
				tgen.T1Gen.push_back(t[i]);
				tgen.T2Gen.push_back(t[N_sample+i]);
			}
		}
		mTGen.push_back(tgen);
		N_tgen=N_sample;
	}


	if(Bound==true)
	{
		double t1=0;
		double t2=0;
		double p1=0;
		double p2=0;
		double con=0;

		for(unsigned int kk=0;kk<mTGen.size();kk++)
		{
			for(unsigned int i=0;i<mTGen.at(kk).T1Gen.size();i++)
			{
				t1=mTGen.at(kk).T1Gen.at(i);
				t2=mTGen.at(kk).T2Gen.at(i);
				if(t1==S_max_x&&t2==S_max_y)
				{
					mTGen.at(kk).T1Gen.at(i)=nan("1");
					mTGen.at(kk).T2Gen.at(i)=nan("2");
				}
				else if(t1==S_max_x)
				{
					p1=opt2dArea(Sdata,0,INFINITY,t2-0.5,t2+0.5,S_max_x,S_max_y);
					p2=opt2dArea(Sdata,t1,INFINITY,t2-0.5,t2+0.5,S_max_x,S_max_y);
					con=-log(p2/p1)/t1;
					rand_val(rand());
					mTGen.at(kk).T1Gen.at(i)=expon(1.0/con)+t1;
				}
				else if(t2==S_max_y)
				{
					p1=opt2dArea(Sdata,t1-0.5,t1+0.5,0,INFINITY,S_max_x,S_max_y);
					p2=opt2dArea(Sdata,t1-0.5,t1+0.5,t2,INFINITY,S_max_x,S_max_y);
					con=-log(p2/p1)/t1;
					rand_val(rand());
					mTGen.at(kk).T2Gen.at(i)=expon(1.0/con)+t2;
				}
			}
		}




	}

}



void GAIT::gait(ExpData* ExpData)
{
	vector<expdata_t> *expdata=ExpData->getExpData();
	vector<double> TGen1_Avg(mTGen.size(),0);
	vector<double> TGen2_Avg(mTGen.size(),0);
	int tempN=0;

	for(unsigned int i=0;i<mTGen.size();i++)
	{
		tempN=N_tgen;
		for(int j=0; j<N_tgen;j++)
		{
			if(mTGen.at(i).T1Gen.at(j)>=0 && mTGen.at(i).T2Gen.at(j)>=0) {TGen1_Avg[i]+=mTGen.at(i).T1Gen.at(j); TGen2_Avg[i]+=mTGen.at(i).T2Gen.at(j);}
			else {tempN=tempN-1; }
		}

		TGen1_Avg[i]=TGen1_Avg[i]/tempN;
		TGen2_Avg[i]=TGen2_Avg[i]/tempN;
	}

	ofstream outfile("pred_T.txt");
	for( unsigned int i=0 ; i<TGen1_Avg.size() ; i++ ) {
		outfile << TGen1_Avg[i]<<"\t"<<TGen2_Avg[i];
		outfile << endl; 
	}
	outfile.close();
}




