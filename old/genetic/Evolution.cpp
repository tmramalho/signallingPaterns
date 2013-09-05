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

Manager *Evolution::train( EvAlgorithm alg , int run_num ) {

	switch (alg) {
		case ALG1:
			return alg_1(run_num);
			break;
		case ALG2:
			return alg_2(run_num);
			break;
		case ALG3:
			return alg_3(run_num);
			break;
		default:
			return NULL;
			break;
	}
}

Manager *Evolution::alg_1( int run_num ) {
		
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
		scores.at(i_agent) =  score/num_to_average;
		if(_verbose) {
			std::cout << "Settler " << i_agent
			<< " ( " << scores.at(i_agent) << " ) ";
			std::cout << std::endl;
		}
	}
	
	//find min and max
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
		
		std::cout << "\nGeneration " << g << std::endl;
		
		if ((g % 10)+1 == 10) {
			genomes_to_file(1,run_num,g+1);
		}
		
		//combination and mutation
		for(int i_agent = 0 ; i_agent < num_agent ; i_agent++) {
			
			next_generation.push_back(new Manager(_pop.at(i_agent)));
			
			if (_do_track > 0) {
				_current_organisms.at(i_agent)->_copy_offspring = new_organism();
				_current_organisms.at(i_agent)->_mutate_offspring = new_organism();
				
				mutation_organisms.push_back( _current_organisms.at(i_agent)->_mutate_offspring );
				_current_organisms.at(i_agent) = _current_organisms.at(i_agent)->_copy_offspring;
			}
			
			if (should_mutate(generator)) {
				if (_do_track > 0) {
					mutation_organisms.at(i_agent)->_did_mutate = true;
				}
				for ( int i_mutation = 0 ; i_mutation < _mutations_per_gen ; i_mutation++ ) {
					next_generation.at(i_agent)->mutate(generator);
				}
			}
			
			double score = 0;
			double score_update = 0;
			for ( int i_run = 0 ; i_run < num_to_average ; i_run++ ) {
				score_update = _ff->run( _pop.at(i_agent) , generator );
				score += score_update;
				if (_do_track > 0) {
					mutation_organisms.at(i_agent)->_evaluations.push_back(score_update);
					if (_do_track == 2) {
						_evaluations.push_back(score_update);
					}
				}
			}
			scores.at(num_agent + i_agent) =  score/num_to_average;
			
		}
		
		// sort
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
	genomes_to_file(1, run_num, -1);
	return _pop.at(i_min);
}

Manager *Evolution::alg_2( int run_num ) {
	
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
		
		if ( (g % 10)+1 == 10) {
			genomes_to_file(2,run_num,g+1);
		}
		
		/* Produce next generation (either track or not) */
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
	genomes_to_file(2, run_num, -1);
	/* Delete Managers we will not return */
	for ( int i_agent = 0 ; i_agent < num_agent ; i_agent++ ) {
		if ( i_agent != i_min ) delete _pop.at(i_agent);
	}
	
	/* Return best */
	return _pop.at(i_min);
	
}

Manager *Evolution::alg_3( int run_num ) {
	
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
		scores.at(i_agent) =  score/num_to_average;
		if(_verbose) {
			std::cout << "Settler " << i_agent
			<< " ( " << scores.at(i_agent) << " ) ";
			std::cout << std::endl;
		}
	}
	
	//find min and max
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
		
		std::cout << "\nGeneration " << g << std::endl;
		
		if ( (g % 10)+1 == 10) {
			genomes_to_file(3,run_num,g+1);
		}
		
		//combination and mutation
		for(int i_agent = 0 ; i_agent < num_agent ; i_agent++) {
			
			next_generation.push_back(new Manager(_pop.at(i_agent)));
			
			if (_do_track > 0) {
				_current_organisms.at(i_agent)->_copy_offspring = new_organism();
				_current_organisms.at(i_agent)->_mutate_offspring = new_organism();
				
				mutation_organisms.push_back( _current_organisms.at(i_agent)->_mutate_offspring );
				_current_organisms.at(i_agent) = _current_organisms.at(i_agent)->_copy_offspring;
			}
			
			if (should_mutate(generator)) {
				if (_do_track > 0) {
					mutation_organisms.at(i_agent)->_did_mutate = true;
				}
				for ( int i_mutation = 0 ; i_mutation < _mutations_per_gen ; i_mutation++ ) {
					next_generation.at(i_agent)->mutate(generator);
				}
			}
			
			double score = 0;
			double score_update = 0;
			for ( int i_run = 0 ; i_run < num_to_average ; i_run++ ) {
				score_update = _ff->run( _pop.at(i_agent) , generator );
				score += score_update;
				if (_do_track > 0) {
					mutation_organisms.at(i_agent)->_evaluations.push_back(score_update);
					if (_do_track == 2) {
						_evaluations.push_back(score_update);
					}
				}
			}
			scores.at(num_agent + i_agent) =  score/num_to_average;
			
		}
		
		// sort
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
	genomes_to_file(3, run_num, -1);
	return _pop.at(i_min);
}


Organism* Evolution::new_organism() {
	Organism* org_ref = new Organism;
	org_ref->_did_mutate = false;
	org_ref->_copy_offspring = NULL;
	org_ref->_mutate_offspring = NULL;
	return org_ref;
}

void Evolution::phylogeny_to_file( std::ofstream& file ) {
	
	int index = -1;
	int generation = 0;
	
	for ( int i_agent = 0 ; i_agent < _na ; i_agent++ ) {
		assign_index( _initial_organisms.at(i_agent) , index );
	}
	
	for ( int i_agent = 0 ; i_agent < _na ; i_agent++ ) {
		phylogeny_to_file( file , _initial_organisms.at(i_agent) , generation );
	}
	
}

void Evolution::evaluations_to_file( std::ofstream& file ) {
	
	int j;
	j = 0;
	
	for ( int i_eval = 0 ; i_eval < _evaluations.size() ; i_eval++ ) {
		file << _evaluations.at(i_eval) << " ";
	}
	
}

vector<Lineage*> Evolution::create_lineages() {
	
	vector<Lineage*> lineages;
	
	for (int i_org = 0 ; i_org < _initial_organisms.size() ; i_org++) {
		lineages.push_back(new Lineage());
		create_lineages(_initial_organisms.at(i_org),lineages,0);
	}
	
	return lineages;
	
}

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

void Evolution::assign_index( Organism* org_ref , int& index) {
	// index is always the index of the last organisms to be added to the genome.
	if (org_ref == NULL) return;
	index++;
	org_ref->_index = index;
	assign_index( org_ref->_copy_offspring , index );
	assign_index( org_ref->_mutate_offspring , index );
	
}

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

void Evolution::delete_organism( Organism* org_ref ) {
	
	if ( org_ref == NULL ) return;
	
	delete_organism( org_ref->_copy_offspring );
	delete_organism( org_ref->_mutate_offspring );
	
	delete org_ref;
	
}

void Evolution::genomes_to_file(int alg_num , int run_num, int gen_num) {
	std::string filename;
	std::ofstream file;
	for ( int i = 0 ; i < _pop.size() ; i++ ) {
		filename = "alg";
		filename += boost::lexical_cast<std::string>(alg_num);
		filename += "/";
		filename += "run";
		filename += boost::lexical_cast<std::string>(run_num);
		filename += "/";
		filename += "genomes";
		filename += boost::lexical_cast<std::string>(gen_num);
		filename += "/genome";
		filename += boost::lexical_cast<std::string>(i);
		filename += ".txt";
		std::cout << filename << "\n";
		file.open(filename.c_str());
		_pop.at(i)->genome_to_file(file,"\t");
		file.close();
	}
}





