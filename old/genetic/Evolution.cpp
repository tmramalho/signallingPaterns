#include "Evolution.h"

Evolution::Evolution(FitnessFunction *ff) {
	SettingsCont *sc_ref = SettingsCont::getInstance();
	_na = sc_ref->_na;
	_ng = sc_ref->_ng;
	_mutations_per_gen = sc_ref->_mutations_per_gen;
	_verbose = sc_ref->_verbose;
	_ff = ff;
	_mutation_probability = sc_ref->_mutation_probability;
	_do_track = sc_ref->_do_track;
	
	if(_verbose) {
		std::cout << "Running differential evolution with"
			<< " n_agents " << _na << std::endl;
	}
}

Evolution::~Evolution() {

	_pop.erase(_pop.begin(),_pop.end());
	
	for ( int i_agent = 0 ; i_agent < _na ; i_agent++ ) {
		delete_organism( _initial_organisms.at(i_agent) );
	}
	
	// The above for loop will also delete the organisms referenced in _current_organisms

}

/* colonize adds Managers to our Manager* vector as the initial population.
 * If we also are in a data tracking mode, it adds the data tracking 
 * objects Organism* to the organisms vectors
 */

void Evolution::colonize( ConstructionMethod m) {
	
	_pop.erase(_pop.begin(),_pop.end());
	
	if (_do_track == 0) {
		for(int i = 0 ; i < _na ; i++) {
			_pop.push_back( new Manager(m) );
		}
	}
	else if (_do_track > 0) {
		for(int i = 0 ; i < _na ; i++) {
			_pop.push_back( new Manager(m) );
			Organism* org_ref = new_organism();
			_initial_organisms.push_back(org_ref);
			_current_organisms.push_back(org_ref);
		}
	}
}

/* train: performs an algorithm based on the method provided to it.
 */

Manager *Evolution::train( EvAlgorithm alg ) {

	switch (alg) {
		case ALG1:
			return alg_1();
			break;
		case ALG2:
			return alg_2();
			break;
		case ALG3:
			return alg_3();
			break;
		default:
			return NULL;
			break;
	}
}

Manager *Evolution::alg_1() {
		
	int num_agent = _pop.size();
	dvec scores(2*num_agent);	
	vector<Manager*> next_generation;
	int num_to_average = 4;
	
	/* Data Tracking */
	vector<Organism*> mutation_organisms;
	
	/* Random Number Generator */
	unsigned seed = std::time(0);
	boost::random::mt19937 generator(seed);
	double probabilities[] = { 1-_mutation_probability , _mutation_probability };
	boost::random::discrete_distribution<> should_mutate(probabilities);

	
	//evaluation
	for(int i_agent = 0 ; i_agent < num_agent ; i_agent++) {
		double score = 0;
		double score_update = 0;
		
		/* Ingrate num_to_average times, and find the score each time */
		for ( int i_run = 0 ; i_run < num_to_average ; i_run++ ) {
			score_update = _ff->run( _pop.at(i_agent) , generator );
			score += score_update;
			if (_do_track > 0) {
				_current_organisms.at(i_agent)->_evaluations.push_back(score_update);
				if (_do_track == 2) {
					_evaluations.push_back(score_update);
				}
			}
		}
		
		/* Average the different integration's scores */
		scores.at(i_agent) =  score/num_to_average;
		
		/* Console output */
		if(_verbose) {
			std::cout << "Settler " << i_agent
			<< " ( " << scores.at(i_agent) << " ) ";
			std::cout << std::endl;
		}
	}
	
	//find min and max in scores vector
	double min = scores.at(0);
	int i_min = 0;
	double max = scores.at(0);
	int i_max = 0;
	for ( int i_agent = 1 ; i_agent < num_agent ; i_agent++ ) {
		if (scores.at(i_agent) < min) {
			min = scores.at(i_agent);
			i_min = i_agent;
		}
		if (scores.at(i_agent) > max) {
			max = scores.at(i_agent);
			i_max = i_agent;
		}
	}
	
	// evolve
	for(int g = 0 ; g < _ng ; g++) {
		
		/* Console ouput */
		std::cout << "\nGeneration " << g << std::endl;
		
		//combination and mutation
		for(int i_agent = 0 ; i_agent < num_agent ; i_agent++) {

			/* Copy parent manager, and store in next_generation vector */
			next_generation.push_back(new Manager(_pop.at(i_agent)));
			
			/* Data Tracking: Create organism representing both the parent genome at this generation, 
			 * and there offspring mutant at this generation
			 */
			if (_do_track > 0) {
				_current_organisms.at(i_agent)->_copy_offspring = new_organism();
				_current_organisms.at(i_agent)->_mutate_offspring = new_organism();
				
				mutation_organisms.push_back( _current_organisms.at(i_agent)->_mutate_offspring );
				_current_organisms.at(i_agent) = _current_organisms.at(i_agent)->_copy_offspring;
			}
			
			/* Decide whetehr to mutate offpring.
			 * If we do mutate, add this information to the offspring mutant Organism.
			 */
			if (should_mutate(generator)) {
				if (_do_track > 0) {
					mutation_organisms.at(i_agent)->_did_mutate = true;
				}
				for ( int i_mutation = 0 ; i_mutation < _mutations_per_gen ; i_mutation++ ) {
					next_generation.at(i_agent)->mutate(generator);
				}
			}
			
			/* For each mutant offspring, calculate the score */
			double score = 0;
			double score_update = 0;
			for ( int i_run = 0 ; i_run < num_to_average ; i_run++ ) {
				score_update = _ff->run( _pop.at(i_agent) , generator );
				score += score_update;
				
				/* Data Tracking: add each evaluation peformed to the mutant offspring's offspring
				 * object.
				 */
				if (_do_track > 0) {
					mutation_organisms.at(i_agent)->_evaluations.push_back(score_update);
					if (_do_track == 2) {
						_evaluations.push_back(score_update);
					}
				}
			}
			/* Average integrations performed. */
			scores.at(num_agent + i_agent) =  score/num_to_average;
			
		}
		
		// sort using scores vector
		// Move both the mutant manager and the offspring organisms to the 
		// _next_generatoin and _current_organisms vector respectively in accordance
		// with the sorting operations.
		for ( int i_new_agent = 0 ; i_new_agent < num_agent ; i_new_agent++ ) {
			if ( scores.at(num_agent+i_new_agent) < max ) {
				int i_replace = i_max;
				Manager *old = _pop.at(i_replace);
				delete old;
				_pop.at(i_replace) = next_generation.at(i_new_agent);
				if (_do_track > 0) {
					_current_organisms.at(i_replace) = mutation_organisms.at(i_new_agent);
				}
				scores.at(i_replace) = scores.at(num_agent+i_new_agent);
				if(_verbose == 2) {
					std::cout << "Changed element " << i_replace
					<< " ( " << scores.at(i_replace) << " ): " << std::endl;
				}
				// Find new max (among first num_agent in scores)
				i_max = 0;
				max = scores.at(i_max);
				for ( int i_agent = 0 ; i_agent < num_agent ; i_agent++ ) {
					if ( scores.at(i_agent) > max ) {
						i_max = i_agent;
						max = scores.at(i_agent);
					}
				}
				// Update min if necessary
				if ( scores.at(i_replace) < min) {
					min = scores.at(i_replace);
					i_min = i_replace;
				}
				std::cout << "Replacement made: max " << max << " min " << min;
				std::cout << " replace " << scores.at(i_replace) << "\n";
				
				
			}
			else { // The i_new_agent mutant doesn't make it into our population
				delete next_generation.at(i_new_agent);
			}
		}
		next_generation.clear();
		mutation_organisms.clear();
		if (min < -9.3) {
			g = _ng;
		}
	}
	return _pop.at(i_min);
}

Manager *Evolution::alg_2( ) {
	
	int num_agent = _pop.size();
	vector<Manager*> next_generation;
	
	/* Data Tracking */
	vector<Organism*> mutation_organisms;
	
	/* Random Number Generator */
	unsigned seed = std::time(0);
	boost::random::mt19937 generator(seed);
	double probabilities[] = { 1-_mutation_probability , _mutation_probability };
	boost::random::discrete_distribution<> should_mutate(probabilities);
	boost::random::uniform_int_distribution<> which_agent(0,num_agent-1);

	
	for ( int g = 0 ; g < _ng ; g++ ) {
		
		std::cout << "\nGeneration " << g << std::endl;
		
		/* Produce next generation (either track or not), and mutate if choose to */
		if (_do_track == 0) {
			for ( int i_agent = 0 ; i_agent < num_agent ; i_agent++ ) {
				next_generation.push_back(new Manager(_pop.at(i_agent)));
				if (should_mutate(generator)) {
					next_generation.at(i_agent)->mutate(generator);//Out here so we only add g to mutation history once
					for ( int i_mutation = 0 ; i_mutation < _mutations_per_gen - 1 ; i_mutation++ ) {
						next_generation.at(i_agent)->mutate(generator);
					}
				}
			}
		}
		/* Data tracking version */
		else if (_do_track > 0) {
			for ( int i_agent = 0 ; i_agent < num_agent ; i_agent++ ) {
				
				_current_organisms.at(i_agent)->_copy_offspring = new_organism();
				_current_organisms.at(i_agent)->_mutate_offspring = new_organism();
				
				mutation_organisms.push_back( _current_organisms.at(i_agent)->_mutate_offspring );
				_current_organisms.at(i_agent) = _current_organisms.at(i_agent)->_copy_offspring;
				
				next_generation.push_back(new Manager(_pop.at(i_agent)));
				if (should_mutate(generator)) {
					mutation_organisms.at(i_agent)->_did_mutate = true;
					for ( int i_mutation = 0 ; i_mutation < _mutations_per_gen ; i_mutation++ ) {
						next_generation.at(i_agent)->mutate(generator);
					}
				}
			}
		}
		
		/* Next generation competes to get into population */
		for ( int i_agent = 0 ; i_agent < num_agent ; i_agent++ ) {
			
			/* Score Calculation */
			int i_target_agent = which_agent(generator); // Agent the old agent will try to displace.
			double new_agent_score = _ff->run(next_generation.at(i_agent), generator);
			double target_agent_score = _ff->run(_pop.at(i_target_agent), generator);
			
			/* Add evaluations to current organisms corresponding to the managers being
			 * evaluated.
			 */
			if ( _do_track > 0 ) {
				_current_organisms.at(i_target_agent)->_evaluations.push_back(target_agent_score);
				mutation_organisms.at(i_agent)->_evaluations.push_back(new_agent_score);
			}
			
			if ( _do_track == 2 ) {
				_evaluations.push_back(new_agent_score);
				_evaluations.push_back(target_agent_score);
			}
			
			/* Replacement of old or death of new */
			if ( new_agent_score < target_agent_score ) {
				
				Manager *old = _pop.at(i_target_agent);
				_pop.at(i_target_agent) = next_generation.at(i_agent);
				delete old;
				
				if (_do_track > 0 ) {
					_current_organisms.at(i_target_agent) = mutation_organisms.at(i_agent);
				}

				std::cout << "Replacement made. New agent " << i_agent << " beat old agent ";
				std::cout << i_target_agent << " " << new_agent_score << " to " << target_agent_score << "\n";
				
				
				/*
				if (new_agent_score < -9.4) {
					new_agent_score = _ff->run(_pop.at(i_target_agent),generator);
					if (new_agent_score < -9.4) {
						g = _ng;
					}
				}
				*/
				/* Return condition */
				if (new_agent_score < -9.4) {
					new_agent_score = _ff->run(_pop.at(i_target_agent), generator);
					if ( _do_track > 0 ) {
						_current_organisms.at(i_target_agent)->_evaluations.push_back(new_agent_score);
					}
					if ( _do_track == 2 ) {
						_evaluations.push_back(new_agent_score);
					}
					if (new_agent_score < -9.4) {
						g = _ng;
					}
				}
				
			}
			else {
				delete next_generation.at(i_agent);
			}
		}
		next_generation.clear();
		mutation_organisms.clear(); // Data tracking... if we aren't tracking data, just clears an already empty vector
	}

	/* Find min among final population */
	int i_min = 0;
	double min = _ff->run(_pop.at(0), generator);
	if ( _do_track > 0 ) {
		_current_organisms.at(0)->_evaluations.push_back(min);
	}
	if ( _do_track == 2 ) {
		_evaluations.push_back(min);
	}
	std::cout << "Final Evaluation: " << i_min << " with score " << min << "\n";
	for ( int i_agent = 1 ; i_agent < num_agent ; i_agent++ ) {
		double next_score = _ff->run(_pop.at(i_agent), generator);
		if ( _do_track > 0 ) {
			_current_organisms.at(i_agent)->_evaluations.push_back(next_score);
		}
		if ( _do_track == 2 ) {
			_evaluations.push_back(next_score);
		}
		if (next_score < min) {
			i_min = i_agent;
			min = next_score;
			std::cout << "Final Evaluation: " << i_min << " with score " << min << "\n";
		}
	}
	/* Delete Managers we will not return */
	for ( int i_agent = 0 ; i_agent < num_agent ; i_agent++ ) {
		if ( i_agent != i_min ) delete _pop.at(i_agent);
	}
	
	/* Return best */
	return _pop.at(i_min);
	
}

Manager *Evolution::alg_3() {
	
	int num_agent = _pop.size();
	dvec scores(2*num_agent);	
	vector<Manager*> next_generation;
	int num_to_average = 2;
	
	/* Data Tracking */
	vector<Organism*> mutation_organisms;
	
	/* Random Number Generator */
	unsigned seed = std::time(0);
	boost::random::mt19937 generator(seed);
	double probabilities[] = { 1-_mutation_probability , _mutation_probability };
	boost::random::discrete_distribution<> should_mutate(probabilities);
	
	
	//evaluation
	for(int i_agent = 0 ; i_agent < num_agent ; i_agent++) {
		double score = 0;
		double score_update = 0;
		
		/* Ingrate num_to_average times, and find the score each time */
		for ( int i_run = 0 ; i_run < num_to_average ; i_run++ ) {
			score_update = _ff->run( _pop.at(i_agent) , generator );
			score += score_update;
			if (_do_track > 0) {
				_current_organisms.at(i_agent)->_evaluations.push_back(score_update);
				if (_do_track == 2) {
					_evaluations.push_back(score_update);
				}
			}
		}
		
		/* Average the different integration's scores */
		scores.at(i_agent) =  score/num_to_average;
		
		/* Console output */
		if(_verbose) {
			std::cout << "Settler " << i_agent
			<< " ( " << scores.at(i_agent) << " ) ";
			std::cout << std::endl;
		}
	}
	
	//find min and max in scores vector
	double min = scores.at(0);
	int i_min = 0;
	double max = scores.at(0);
	int i_max = 0;
	for ( int i_agent = 1 ; i_agent < num_agent ; i_agent++ ) {
		if (scores.at(i_agent) < min) {
			min = scores.at(i_agent);
			i_min = i_agent;
		}
		if (scores.at(i_agent) > max) {
			max = scores.at(i_agent);
			i_max = i_agent;
		}
	}
	
	// evolve
	for(int g = 0 ; g < _ng ; g++) {
		
		/* Console ouput */
		std::cout << "\nGeneration " << g << std::endl;
		
		//combination and mutation
		for(int i_agent = 0 ; i_agent < num_agent ; i_agent++) {
			
			/* Copy parent manager, and store in next_generation vector */
			next_generation.push_back(new Manager(_pop.at(i_agent)));
			
			/* Data Tracking: Create organism representing both the parent genome at this generation, 
			 * and there offspring mutant at this generation
			 */
			if (_do_track > 0) {
				_current_organisms.at(i_agent)->_copy_offspring = new_organism();
				_current_organisms.at(i_agent)->_mutate_offspring = new_organism();
				
				mutation_organisms.push_back( _current_organisms.at(i_agent)->_mutate_offspring );
				_current_organisms.at(i_agent) = _current_organisms.at(i_agent)->_copy_offspring;
			}
			
			/* Decide whetehr to mutate offpring.
			 * If we do mutate, add this information to the offspring mutant Organism.
			 */
			if (should_mutate(generator)) {
				if (_do_track > 0) {
					mutation_organisms.at(i_agent)->_did_mutate = true;
				}
				for ( int i_mutation = 0 ; i_mutation < _mutations_per_gen ; i_mutation++ ) {
					next_generation.at(i_agent)->mutate(generator);
				}
			}
			
			/* For each mutant offspring, calculate the score */
			double score = 0;
			double score_update = 0;
			for ( int i_run = 0 ; i_run < num_to_average ; i_run++ ) {
				score_update = _ff->run( _pop.at(i_agent) , generator );
				score += score_update;
				
				/* Data Tracking: add each evaluation peformed to the mutant offspring's offspring
				 * object.
				 */
				if (_do_track > 0) {
					mutation_organisms.at(i_agent)->_evaluations.push_back(score_update);
					if (_do_track == 2) {
						_evaluations.push_back(score_update);
					}
				}
			}
			/* Average integrations performed. */
			scores.at(num_agent + i_agent) =  score/num_to_average;
			
		}
		
		// sort using scores vector
		// Move both the mutant manager and the offspring organisms to the 
		// _next_generatoin and _current_organisms vector respectively in accordance
		// with the sorting operations.
		for ( int i_new_agent = 0 ; i_new_agent < num_agent ; i_new_agent++ ) {
			if ( scores.at(num_agent+i_new_agent) < max ) {
				int i_replace = i_max;
				Manager *old = _pop.at(i_replace);
				delete old;
				_pop.at(i_replace) = next_generation.at(i_new_agent);
				if (_do_track > 0) {
					_current_organisms.at(i_replace) = mutation_organisms.at(i_new_agent);
				}
				scores.at(i_replace) = scores.at(num_agent+i_new_agent);
				if(_verbose == 2) {
					std::cout << "Changed element " << i_replace
					<< " ( " << scores.at(i_replace) << " ): " << std::endl;
				}
				// Find new max (among first num_agent in scores)
				i_max = 0;
				max = scores.at(i_max);
				for ( int i_agent = 0 ; i_agent < num_agent ; i_agent++ ) {
					if ( scores.at(i_agent) > max ) {
						i_max = i_agent;
						max = scores.at(i_agent);
					}
				}
				// Update min if necessary
				if ( scores.at(i_replace) < min) {
					min = scores.at(i_replace);
					i_min = i_replace;
				}
				std::cout << "Replacement made: max " << max << " min " << min;
				std::cout << " replace " << scores.at(i_replace) << "\n";
				
				
			}
			else { // The i_new_agent mutant doesn't make it into our population
				delete next_generation.at(i_new_agent);
			}
		}
		next_generation.clear();
		mutation_organisms.clear();
		if (min < -9.3) {
			g = _ng;
		}
	}
	return _pop.at(i_min);
}

/* new_organism()
 * --------------------------------------------------------------------------------------------
 * Useful only for data tracking. Creates a new organism object on the heap and initialize
 * its data members as follows:
 *
 *		did_mutate = false because the organism hasn't been mutate at initiation
 *		_copy_offspring, _mutate_offspring = NULL because it has no offspring yet
 */
Organism* Evolution::new_organism() {
	Organism* org_ref = new Organism;
	org_ref->_did_mutate = false;
	org_ref->_copy_offspring = NULL;
	org_ref->_mutate_offspring = NULL;
	return org_ref;
}
						 
/* phylogeny_to_file(file)
 * --------------------------------------------------------------------------------------------
 * File Format: List of
 *
 *					000 <Organism i Index>
 *					001 <Organism i Generation>
 *					002	<Organism i did_mutate>
 *					003 <Eval 1> < Eval 2 > ... < Eval n>
 *					004 <Index of Copy Offspring> 
 *					005 <Index of Mutate Offspring
 *					
 * First it must assign a unique index to each organism in the phylogeny. It does this by
 * calling the assign_index method on each organism in the _initial_organisms vector, starting
 * the index count at -1 (see method for info on implementation.
 *
 * Second, it outputs each Organism's info to the file, in succession, by calling the
 * phylogeny_to_file method on each Organism in the _initial_organisms vector. This method
 * outputs to file an organism and all its descedents, in an order determined by the recursive 
 * path taken in traversing the phylogeny tree, a path which will match the order index are 
 * assigned by the assign_index method, which also works through recursive traversal.
 */
void Evolution::phylogeny_to_file( std::ofstream& file ) {
	
	int index = -1;
	int generation = 0;
	
	/* First give each organism an index */
	for ( int i_agent = 0 ; i_agent < _na ; i_agent++ ) {
		assign_index( _initial_organisms.at(i_agent) , index );
	}
	
	/* Now sequentially output each organism's info the file */
	/* We loop through the agents in the first generation. Because there is not combination,
	 * lineages never combine, so we call the phylogeny_to_file function that prints an organism
	 * to the file, and all the organisms deriving from it at some later point in the evolutionary
	 * history.
	 */
	for ( int i_agent = 0 ; i_agent < _na ; i_agent++ ) {
		phylogeny_to_file( file , _initial_organisms.at(i_agent) , generation );
	}
	
}
						 
/* Outputs list of all evaluations to file */

void Evolution::evaluations_to_file( std::ofstream& file ) {
	
	int j;
	j = 0;
	
	for ( int i_eval = 0 ; i_eval < _evaluations.size() ; i_eval++ ) {
		file << _evaluations.at(i_eval) << " ";
	}
	
}

/* create_lineages
 * --------------------------------------------------------------------------------------------
 * Returns a vector pointed to all the lineages in our entire phylogeny. The entire phylogeny
 * is the tree consisting of each initial organism and all their offpsring.
 *
 * Achieves this by calling the recursive create_lineages method on each initial organism. This
 * recursive method adds to the vector passed to it by reference all lineages in whcih org_ref
 * takes part.
 *
 * See create_lineages(org_ref,lineages,generation) below for more detail.
 */

vector<Lineage*> Evolution::create_lineages() {
	
	vector<Lineage*> lineages;
	
	for (int i_org = 0 ; i_org < _initial_organisms.size() ; i_org++) {
		lineages.push_back(new Lineage());
		create_lineages(_initial_organisms.at(i_org),lineages,0);
	}
	
	return lineages;
	
}
						 
/* phylogeny_to_file(file,org_ref,gen)
 * --------------------------------------------------------------------------------------------
 * Recursive function.
 *
 * Input Conditions: 
 *		a) The Organism pointed to by org_ref has not been output to the file yet.
 *		b) gen is the generation of the Organism pointed to by org_ref
 *		c) org_ref may be NULL
 *
 * Method Function:
 *		Outputs the organism passed to it to the file, as well as all of the organisms
 *		descendants. It achieves this by calling itself on each of the organisms offspring.
 */

void Evolution::phylogeny_to_file( std::ofstream& file , Organism* org_ref , int gen ) {
	
	if (org_ref == NULL) {
		return;
	}
	
	file << "000 " << org_ref->_index << "\n";
	file << "001 " << gen << "\n";
	file << "002 " << org_ref->_did_mutate << "\n";
	file << "003 ";
	for ( int i_eval = 0 ; i_eval < org_ref->_evaluations.size() ; i_eval++ ) {
		file << org_ref->_evaluations.at(i_eval) << " ";
	}
	file << "\n";
	
	if (org_ref->_copy_offspring == NULL) {
		file << "004 " << -1 << "\n";
		file << "005 " << -1 << "\n" << "\n";
	}
	else {
		file << "004 " << org_ref->_copy_offspring->_index << "\n";
		file << "005 " << org_ref->_mutate_offspring->_index << "\n" << "\n";
	}
	phylogeny_to_file( file , org_ref->_copy_offspring , gen+1 );
	phylogeny_to_file( file , org_ref->_mutate_offspring , gen+1 );
	
}

/* assign_indes(org_ref,index)
 * --------------------------------------------------------------------------------------------
 * Recursive function
 *
 * Input Conditions:
 *		a) org_ref points to an organism that has not yet been assigned an index.
 *		b) index, passed by reference, is the index of the LAST assigned index. Thus, it must
 *		be incremented before being assigned to the Organism pointed to by org_ref.
 *		c) org_ref may be NULL
 *
 * Method Function:
 *		Assigns to the Organism pointed to by org_ref the index index+1, and assigns an index
 *		unique index to each of *org_ref's offspring. The indices assigned will form a set
 *		of consecutive integers (though the order in which they are assigned depends on the 
 *		recursive traversal path taken... described below). Thus, when the function reterns,
 *		index, passed by reference has been incremented to the last index assigned in the 
 *		entire subtree of the phylogeny graph descending from *org_ref.
 * 
 * The order of the recursive traversal leads to indices assigned to copy offspring first, 
 * represented by the left-most branches in the example tree below, which would result if the 
 * function had been called at the top node with index = 6.
 *
 *											7
 *										   / \
 *										  8  15
 *                                       / \
 *                                      9  10
 *										  /  \
 *										 11  14
 *										/  \ 
 *									   12  13
 *                                         
 * When the funtion returns, in this example case, index = 15.
 */
void Evolution::assign_index( Organism* org_ref , int& index) {
	// index is always the index of the last organisms to be added to the genome.
	if (org_ref == NULL) return;
	index++;
	org_ref->_index = index;
	assign_index( org_ref->_copy_offspring , index );
	assign_index( org_ref->_mutate_offspring , index );
	
}

						 
/* create_lineages(org_ref,lineages,generation
 * --------------------------------------------------------------------------------------------
 * Recursive function.
 *
 * Input Conditions:
 *		a) The last element of the lineages vector points to a lineage going from the root 
 *		Organism to the parent Organism of org_ref. The information about org_ref (ie. any 
 *		evaluations performed on it or whether or not it has mutated, has not yet been added
 *		to the lineage. Further, this means lineaves CANNOT be an empty vector.
 *		b) org_ref is NOT NULL
 *		c) generation is the generation of *org_ref.
 *
 * Method Function
 *		When the function returns, the vector lineages, passed by reference, will contain all
 *		the lineages that descend from org_ref. The last lineage in the vector will be a 
 *		complete lineage that passes the org_ref, in particular, the one that was added last.
 *
 * Recursion
 *		1. First it adds to the last lineage in the lineages vector the evaluations that
 *		occured to this Organism, and whether or not it mutated. This works because the last 
 *		lineage in the vector does not yet have this info but is the lineage that leads to
 *		Organism's parent (see a).
 *		2. This is not the end of a lineage only if org_ref's offspring are not NULL. Thus, we
 *		check this, and if they are NULL, there is nothing more to do.
 *		3. We copy the current lineage. All lineages that go through org_ref will be the same
 *		in all the data up to this point.
 *		4. Call create_lineages on the copy_offspring and the current lineages vector. If the
 *		function operates correctly, this will add to the lineages vector all the lineages that
 *		pass from org_ref to copy_offspring.
 *		5. When 4 is over, the last lineage in the lineages vector is the most recently added
 *		lineage that passes through copy_offspring. We want to add all lineages that pass
 *		through mutate_offspring now. But, to maintain the input conditions, we must call
 *		create_lineages on the mutate_offspring only when the lineages vector's last element
 *		is a lineage ending at org_ref (the parent of mutate_offspring). We add the copy
 *		we made in 3 so that this is the case.
 */
void Evolution::create_lineages( Organism* org_ref , std::vector<Lineage*>& lineages , int generation ) {
	
	// We assume org_ref is not NULL
	
	int i_curr_lineage = lineages.size()-1;

	for (int i_eval = 0; i_eval < org_ref->_evaluations.size(); i_eval++) {
		lineages.at(i_curr_lineage)->add_evaluation(org_ref->_evaluations.at(i_eval),generation);
	}
	if (org_ref->_did_mutate == true) {
		lineages.at(i_curr_lineage)->add_mutation(generation);
	}
	
	if (org_ref->_copy_offspring != NULL) {
		Lineage *lin_copy = new Lineage(lineages.at(lineages.size()-1));
		create_lineages(org_ref->_copy_offspring,lineages,generation+1);
		lineages.push_back(lin_copy); // new lineage for the new path
		create_lineages(org_ref->_mutate_offspring,lineages,generation+1);
	}
	
}

/* delete_organism
 * --------------------------------------------------------------------------------------------
 * Recursive Function
 *
 * Input Conditions:
 *		a) org_ref points to an organism that has not yet been deleted, or is NULL
 *
 * Method Function:
 *		Delete all of *org_ref's descendents as well as org_ref itself.
 *
 * Recursive Call: 
 *		Call delete_organism on each descendand, and then delete org_ref itself.
 */
void Evolution::delete_organism( Organism* org_ref ) {
	
	if ( org_ref == NULL ) return;
	
	delete_organism( org_ref->_copy_offspring );
	delete_organism( org_ref->_mutate_offspring );
	
	delete org_ref;
	
}




