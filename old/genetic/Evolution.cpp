#include "Evolution.h"

Evolution::Evolution(FitnessFunction *ff) {
	SettingsCont *sc_ref = SettingsCont::getInstance();
	_na = sc_ref->_na;
	_ng = sc_ref->_ng;
	_mutations_per_gen = sc_ref->_mutations_per_gen;
	_verbose = sc_ref->_verbose;
	_ff = ff;
	
	if(_verbose) {
		std::cout << "Running differential evolution with"
			<< " n_agents " << _na << std::endl;
	}
}

Evolution::~Evolution() {

	_pop.erase(_pop.begin(),_pop.end());

}

void Evolution::colonize( ConstructionMethod m) {
	
	_pop.erase(_pop.begin(),_pop.end());
	
	for(int i = 0 ; i < _na ; i++) {
		_pop.push_back( new Manager(m) );
	}
	
}

Manager *Evolution::train() {

	int size = _pop.size();
	dvec *scores = new dvec(size);
	Manager *child;
	
	/* Random Number Generator */
	unsigned seed = std::time(0);
	boost::random::mt19937 generator(seed);

	//evaluation
	for(int i = 0 ; i < size ; i++) {
		scores->at(i) =  _ff->run( _pop.at(i) , generator );
		if(_verbose) {
			std::cout << "Settler " << i
				<< " ( " << scores->at(i) << " ) ";
			std::cout << std::endl;
		}
	}
	
	double min = scores->at(0);
	int i_min = 0;
	double max = scores->at(0);
	int i_max = 0;
	
	for ( int i = 1 ; i < size ; i++ ) {
		if (scores->at(i) < min) {
			min = scores->at(i);
			i_min = i;
		}
		if (scores->at(i) > max) {
			max = scores->at(i);
			i_max = i;
		}
	}

	for(int g = 0 ; g < _ng ; g++) {
		
		std::cout << "Generation " << g << std::endl;

		//combination and mutation
		for(int i = 0 ; i < size ; i++) {
			
			child = new Manager(_pop.at(i));
			
			for ( int j = 0 ; j < _mutations_per_gen ; j++ ) {
				child->mutate(generator);
				if (child->has_comb_reac()) {
					int l = 2;
				}
			}
			
			double child_score = _ff->run(child,generator);

			if( child_score < max) {
				
				Manager *old = _pop.at(i_max);
				_pop.at(i_max) = child;
				scores->at(i_max) = child_score;
				delete old;
	
				if(_verbose == 2) {
					std::cout << "Changed element " << i
						<< " ( " << scores->at(i) << " ): " << std::endl;
				}
				
				if(child_score < min) {
					min = child_score;
					i_min = i_max; // Then child_score has been put in this position (child_score < max too).
					if(_verbose) std::cout << "New best: " << min << std::endl;
				}
				
				/* Find new max */
				i_max = 0;
				max = scores->at(i_max);
				for ( int i = 0 ; i < scores->size() ; i++ ) {
					if ( scores->at(i) > scores->at(i_max) ) {
						i_max = i;
						max = scores->at(i);
					}
				}
				
			}
		}
	}

	delete scores;

	return _pop.at(i_min);
}

