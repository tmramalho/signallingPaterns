#ifndef TESTFUNCTION_H_
#define TESTFUNCTION_H_

#include "FitnessFunction.h"
#include <math.h>

class TestFunction : public FitnessFunction
{
public:
	TestFunction(int dim, int type);
	~TestFunction();
	
	double run(dvec *x);
	
private:
	int _dim, _type;
};

#endif /*TESTFUNCTION_H_*/
