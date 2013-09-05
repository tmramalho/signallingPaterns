/*
 *  TempFF.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 8/12/13.
 *
 */


#ifndef TEMPFF_H
#define TEMPFF_H

# include "FitnessFunction.h"

class TempFF : public FitnessFunction
{
public:
	TempFF() {}
	virtual ~TempFF() {}
	
	virtual double run(Manager *x,boost::random::mt19937& generator);
	
};

#endif