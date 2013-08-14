#ifndef FITNESSFUNCTION_H_
#define FITNESSFUNCTION_H_

#include <vector>
#include <string>
#include "../numeric/dvec.h"
#include "../../Manager.h"

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
	
	virtual double run(Manager *x)=0;
};

#endif /*FITNESSFUNCTION_H_*/
