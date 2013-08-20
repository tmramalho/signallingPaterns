#ifndef _EVOLUTION_H
#define _EVOLUTION_H

#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
#include "../numeric/dvec.h"
#include "../opt/FitnessFunction.h"
#include "../opt/DNScore.h"
#include "../helpers/SettingsCont.h"
#include "../../Manager.h"
#include "../../ConstructionMethod.h"
#include "omp.h"

# include <boost/random/uniform_int_distribution.hpp>
# include <boost/random/uniform_real_distribution.hpp>
# include <boost/random/mersenne_twister.hpp>

class Evolution {

public:
	Evolution(FitnessFunction *ff);

	virtual ~Evolution();

	void colonize( ConstructionMethod m );

	Manager *train();
		
private:
	FitnessFunction *_ff;
	std::vector<Manager*> _pop;
	int _verbose;
	int _ng;
	int _na;
	int _mutations_per_gen;
};

#endif
