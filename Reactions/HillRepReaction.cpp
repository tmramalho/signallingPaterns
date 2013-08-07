/*
 *  HillDegReaction.cpp
 *  Signalling Patterns
 *
 *  Created by Michael Celentano 2 on 8/5/13.
 *
 */

#include "HillRepReaction.h"

HillRepReaction::HillRepReaction() {
	
	_type = HILL_REPRESSION;
	
	_num_part = 1;
	
}

HillRepReaction::HillRepReaction( int i_repressor , int i_repressed ,
								   double kinetic , double K , double cooperativity ) {
	
	_type = HILL_REPRESSION;
	
	_num_part = 1;
	
	_i_repressor = i_repressor;
	_i_repressed = i_repressed;
	_kinetic = kinetic;
	_K = K;
	_cooperativity = cooperativity;
	
}

HillRepReaction::HillRepReaction(HillRepReaction* newOne) {
	
	_type = HILL_REPRESSION;
	
	_num_part = 1;
	
	_i_repressor = newOne->_i_repressor;
	_i_repressed = newOne->_i_repressed;
	_kinetic = newOne->_kinetic;
	_K = newOne->_K;
	_cooperativity = newOne->_cooperativity;
	
}

HillRepReaction::~HillRepReaction() {}

int HillRepReaction::get_i_part( int part_num ) {
	switch (part_num) {
		case 0:
			return _i_repressed;
			break;
		default:
			return NEXIST;
			break;
	}
}

void HillRepReaction::react( dmat& curr_tissue , dmat& dx_dt ,
							 std::vector< std::vector<int>* >& neighbors ,
							 int i_curr_cell ) {
	
	double r = 
	_kinetic
	* _K
	/
	( _K + pow( curr_tissue.at(i_curr_cell , _i_repressor ) , _cooperativity ));
	
	dx_dt.at(i_curr_cell,_i_repressed) += r;
}

void HillRepReaction::update_indices( int first_index , int num_insertion ) {
	if (_i_repressor >= first_index) _i_repressor += num_insertion;
	if (_i_repressed >= first_index) _i_repressed += num_insertion;
}

void HillRepReaction::mutate ( boost::random::mt19937& generator ) {
	
	boost::random::uniform_int_distribution<> which_kinetic(0,2);
	boost::random::uniform_real_distribution<> mutation_factor(0.0,2.0);
	
	int kinetic_num = which_kinetic(generator);
	double mut_factor = mutation_factor(generator);
	
	switch ( kinetic_num ) {
		case 0:
			_kinetic *= mut_factor;
			break;
		case 1:
			_K *= mut_factor;
			break;
		case 2:
			_cooperativity *= mut_factor;
			break;
		default:
			break;
	}
	
}

void HillRepReaction::print_info ( std::string line_start ) {
	
	std::cout << line_start << "ReactionType: Hill Repression Reaction" << std::endl;
	//std::cout << line_start << "_num_part: " << _num_part << std::endl;
	std::cout << line_start << "Index of Repressor Protein: " << _i_repressor << std::endl;
	std::cout << line_start << "Index of Repressed Protein: " << _i_repressed << std::endl;
	std::cout << line_start << "Kinetic: " << _kinetic << std::endl;
	std::cout << line_start << "K: " << _K << std::endl;
	std::cout << line_start << "Cooperativity: " << _cooperativity << std::endl;
	
}




