

#ifndef _random_h
#define _random_h

//# random number generators are defined in this file

#include <random>
#include <cmath>

#define inv_sqrt_2pi 0.3989422804014327

//////poisson prob den fun
long double dpoi(const long double k, const long double lambda) {
	return exp(k * log(lambda) - lgamma(k + 1.0) - lambda);
}

////gaussian prob den fun
long double dgaus(const long double x, const long double m, const long double s){
    double a = (x - m) / s;
    return inv_sqrt_2pi / s * exp(-0.5f * a * a);
}

////binomial prob d f
double PDF(double p, int s, int n){
  // assert(p>=0 && p<=1);
  double m=1;
  for(int i=0; i<n-s; i++)
    m*=(s+i+0.0)/i;
  m*=pow(p,s)*pow(1-p,n-s);
  return m;
}

using namespace std;

//#generall class for several number generators based on GNU libs
class NumGenerator {

	default_random_engine gen; //standard number generator
	binomial_distribution<int> d; //binomial generator
	uniform_real_distribution<float> u; //uniform generator
	normal_distribution<float> nd; //normal generator
	gamma_distribution<float> ga; //gamma generator
	cauchy_distribution<float> cau; //cauchy generator

	public:

		NumGenerator(long int seed){
			gen= default_random_engine(seed);
		}

		int binGen(int iBE, float rate){
			// iBE - number of indep Bernoulli-distributed experiments each generated value is said to simulate
			d= binomial_distribution<int>(iBE,rate);
			return d(gen);
		}

		float uniGen(){ //[a, b), a=0, b=1
			return u(gen);
		}

		float gammGen(float mean, float theta){
			//note: the docs on this c++ function are confusing and define gamm dist as using alpha and beta
			//however, by the resulting mean, it is clear that this function is actually the version of the
			//gamm distribution which is implemented with parameters k and theta, with mean= k*theta
			ga= gamma_distribution<float>(mean/theta, theta);
			return ga(gen);
		}

		float norGen(float mean, float stdev){
			nd= normal_distribution<float>(mean, stdev);
			return nd(gen);
		}

		float cauGen(float peak, float scale){
			cau= cauchy_distribution<float>(peak, scale);
			return cau(gen);
		}

} ;

NumGenerator* numGen;


#endif
