#ifndef _DIFFEVOLUTION_H
#define _DIFFEVOLUTION_H

#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
#include "../numeric/dvec.h"
#include "../opt/FitnessFunction.h"
#include "../helpers/SettingsCont.h"
#include "mtwist.h"
#include "omp.h"

class DiffEvolution {

	public:
		DiffEvolution(FitnessFunction *ff);

		virtual ~DiffEvolution();

		void copyAgent(int pos, dvec *x);
		void colonize(dvec *x);

		dvec *train();
		
	private:
		FitnessFunction *_ff;
		std::vector<dvec*> *_pop;
		std::vector<dvec*> *_next;
		int _verbose;
		int _acc;
		int _ng;
};

#endif
