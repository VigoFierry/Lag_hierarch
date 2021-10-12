// Generating random variables
// Date     : March 19th 2017

/* #include <iostream>
#include <random>  // package for generating random variables - http://en.cppreference.com/w/cpp/numeric/random or http://itscompiling.eu/2016/04/11/generating-random-numbers-cpp/
*/

#include "Header.h"



	
double rnd()  // This function returns a random floating point number between 0 and 1 - NEPOUZIVAT RADSI
{ 
	srand(time(NULL));
	return double(rand()) / RAND_MAX; 
} 

double uniform(int a, int b)
{
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

	std::uniform_real_distribution<> dis(a, b);
	return dis(gen);

}

int uniform_int(int a, int b)
{
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

	std::uniform_int_distribution<> dis(a, b);
	return dis(gen);

}

double normal(double mu, double sigma) // 
{
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

	std::normal_distribution<> d(mu, sigma);
	return d(gen);
}

double poisson(int lambda)
{
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

	std::poisson_distribution<> d(lambda);
	return d(gen);

}

double gamma(double alfa, double beta)
{
	double rate = 1 / beta;
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()

	std::gamma_distribution<> d(alfa,rate);
	return d(gen);

}

double triangle(double a, double b, double c)		// triangle distribution
{
	double F = (c - a) / (b - a);
	double r = uniform(0, 1);						// using inverse tranformation method
	if ((a < c) & (c < b)) { 
		if (r < F) { return (a + sqrt(r*(b - a)*(c - a))); }
		else { return (b - sqrt((1 - r)*(b - a)*(b - c))); }
	}
	else { std::cout << "Please enter correct parameters (i.e. a < c < b) \n"; return -1; }
}

double exponential(double lambda)
{
	double r = uniform(0, 1);						// using inverse tranformation method

	return (- log(r)/lambda);
}



// uniform distribution on sphere, cf. http://corysimon.github.io/articles/uniformdistn-on-sphere/
void sphere_uniform(double &s1, double &s2, double &s3) {

	// We have to be careful in the case that the vector has a norm close to zero, in which we must worry about floating point precision by dividing by a very small number.This is the reason for the while loop.
	
	s1 = 0; s2 = 0; s3 = 0;  // initialize so we go into the while loop

	double norm = 0;

	while (norm < .0001) {
		s1 = normal(0, 1);  // random standard normal
		s2 = normal(0, 1);
		s3 = normal(0, 1);

		norm = sqrt(pow(s1, 2) + pow(s2, 2) + pow(s3, 2));
	}

	s1 = s1 / norm;  // normalize to unit norm
	s2 = s2 / norm;
	s3 = s3 / norm;

}