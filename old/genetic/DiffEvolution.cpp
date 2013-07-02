#include "DiffEvolution.h"

DiffEvolution::DiffEvolution(FitnessFunction *ff) {
	SettingsCont *s = SettingsCont::getInstance();
	int nPoints = s->dim;
	int nAgents = s->na;
	_ng = s->ng;
	_acc = s->acc;
	_verbose = s->verbose;
	_ff = ff;
	_pop = new std::vector<dvec*>;
	_next = new std::vector<dvec*>;
	
	if(_verbose) {
		std::cout << "Running differential evolution with"
			<< " nAgents " << nAgents << std::endl;
	}

	for(int i = 0 ; i < nAgents ; i++) {
		dvec *agent = new dvec(nPoints);
		for(int j = 0 ; j < nPoints ; j++) {
			agent->at(j) = mt_drand();
		}
		_pop->push_back(agent);
		_next->push_back(new dvec(nPoints));
	}

	omp_set_num_threads(s->numThreads);

}
DiffEvolution::~DiffEvolution() {
	for(uint i = 0 ; i < _next->size() ; i++) {
		delete _next->at(i);
	}

	delete _next;
	for(uint i = 0 ; i < _pop->size() ; i++) {
		delete _pop->at(i);
	}
		
	delete _pop;
}

void DiffEvolution::copyAgent(int pos, dvec *x) {
	uint size = _pop->at(pos)->size();
	if(x->size() >= size) {
		for (uint i = 0; i < size; i++) {
			_pop->at(pos)->at(i) = x->at(i);
		}
	}
}

void DiffEvolution::colonize(dvec *x) {
	int nAgents = _pop->size();
	int nPoints = x->size();
	for(int i = 0 ; i < nAgents ; i++) {
		dvec *agent = _pop->at(i);
		for(int j = 0 ; j < nPoints ; j++) {
			if(i == 0) agent->at(j) = x->at(j);
			else agent->at(j) = x->at(j) + mt_drand() - 0.5;
		}
	}
}

dvec *DiffEvolution::train() {

	uint size = _pop->size();
	uint dsize = _pop->at(0)->size();
	dvec *scores = new dvec(size);

	double min = 1e100;
	double minPos = 0;
	double cr = 0.9;
	double f = 0.8;

	if(size < 4) return NULL;

	//evaluation
	for(uint i = 0 ; i < size ; i++) {
		scores->at(i) =  _ff->run(_pop->at(i));
		if(_verbose) {
			std::cout << "Settler " << i
				<< " ( " << scores->at(i) << " ) ";
			if(_verbose == 2) std::cout << *(_pop->at(i));
			std::cout << std::endl;
		}
	}

	for(int g = 0 ; g < _ng ; g++) {

		std::cout << "Generation " << g << std::endl;

		//combination and mutation
#pragma omp parallel for
		for(uint i = 0 ; i < size ; i++) {
			uint jmin = mt_lrand() % dsize;
			uint r1 = mt_lrand() % size;
			uint r2 = mt_lrand() % size;
			uint r3 = mt_lrand() % size;
			while (r1 == i) r1 = mt_lrand() % size;
			while (r2 == r1 || r2 == i) r2 = mt_lrand() % size;
			while (r3 == r1 || r3 == r1 || r3 == i) r3 = mt_lrand() % size;
			for(uint j = 0 ; j < dsize ; j++) {
				if (j == jmin || mt_drand() < cr) {
					double ng = _pop->at(r1)->at(j) +
										f * (_pop->at(r2)->at(j)-_pop->at(r3)->at(j));
					_next->at(i)->at(j) = ng;
				} else {
					_next->at(i)->at(j) = _pop->at(i)->at(j);
				}
			}
			
			double childScore = _ff->run(_next->at(i));

#pragma omp critical
			{
				if( childScore < scores->at(i)) {
					_pop->at(i)->copy(*(_next->at(i)));
					scores->at(i) = childScore;
					if(_verbose == 2) {
						std::cout << "Changed element " << i
							<< " ( " << scores->at(i) << " ): "
							<< *(_pop->at(i)) << std::endl;
					}
					if(childScore < min) {
						min = childScore;
						minPos = i;
						if(_verbose) std::cout << "New best: " << min << std::endl;
					}
				}
			}
		}
	}

	delete scores;
	
	dvec *result = new dvec(dsize);
	result->copy(*(_pop->at(minPos)));

	return result;
}

