#ifndef FITNESSFUNCTION_H_
#define FITNESSFUNCTION_H_

#include <vector>
#include <string>
#include "../numeric/dvec.h"
#include "../../Manager.h"

/* Boost Random Libraries */
# include <boost/random/uniform_int_distribution.hpp>
# include <boost/random/uniform_real_distribution.hpp>
# include <boost/random/mersenne_twister.hpp>
# include <boost/random/normal_distribution.hpp>

#define PI 3.14159265
#define NOISE 0.00155231
#define NOISE2 0.01
#define SIGMA0 0.00165431
#define PVAL 2.74861

class FitnessFunction
{
public:
	FitnessFunction() {}
	virtual ~FitnessFunction() {}
	
	virtual double run(Manager *x,boost::random::mt19937& generator)=0;
};

#endif /*FITNESSFUNCTION_H_*/
