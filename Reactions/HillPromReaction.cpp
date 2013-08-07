/*
 *  HillPromReaction.cpp
 *	Signalling Patterns
 *
 *  Created by Michael Celentano 2 on 8/5/13.
 *
 */

#include "HillPromReaction.h"

HillPromReaction::HillPromReaction() {
	
	_type = HILL_PROMOTION;
	
	_num_part = 1;

}

HillPromReaction::HillPromReaction( int i_promoter , int i_promoted ,
								   double kinetic , double K , double cooperativity ) {
	
	_type = HILL_PROMOTION;
	
	_num_part = 1;
	
	_i_promoter = i_promoter;
	_i_promoted = i_promoted;
	_kinetic = kinetic;
	_K = K;
	_cooperativity = cooperativity;
	
}

HillPromReaction::HillPromReaction(HillPromReaction* newOne) {
	
	_type = HILL_PROMOTION;
	
	_num_part = 1;
	
	_i_promoter = newOne->_i_promoter;
	_i_promoted = newOne->_i_promoted;
	_kinetic = newOne->_kinetic;
	_K = newOne->_K;
	_cooperativity = newOne->_cooperativity;
	
}

HillPromReaction::~HillPromReaction() {}

int HillPromReaction::get_i_part( int part_num ) {
	switch (part_num) {
		case 0:
			return _i_promoted;
			break;
		default:
			return NEXIST;
			break;
	}
}

void HillPromReaction::react( dmat& curr_tissue , dmat& dx_dt ,
							 std::vector< std::vector<int>* >& neighbors ,
							 int i_curr_cell ) {

	double r = 
	_kinetic
	* pow( curr_tissue.at(i_curr_cell , _i_promoter ) , _cooperativity )
	/
	( _K + pow( curr_tissue.at(i_curr_cell , _i_promoter ) , _cooperativity ));
	 
	dx_dt.at(i_curr_cell,_i_promoted) += r;
}

void HillPromReaction::update_indices( int first_index , int num_insertion ) {
	if (_i_promoter >= first_index) _i_promoter += num_insertion;
	if (_i_promoted >= first_index) _i_promoted += num_insertion;
}

void HillPromReaction::mutate ( boost::random::mt19937& generator ) {

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

void HillPromReaction::print_info ( std::string line_start ) {
	
	std::cout << line_start << "ReactionType: Hill Promotion Reaction" << std::endl;
	//std::cout << line_start << "_num_part: " << _num_part << std::endl;
	std::cout << line_start << "Index of Promoting Protein: " << _i_promoter << std::endl;
	std::cout << line_start << "Index of Promoted Protein: " << _i_promoted << std::endl;
	std::cout << line_start << "Kinetic Rate: " << _kinetic << std::endl;
	std::cout << line_start << "K: " << _K << std::endl;
	std::cout << line_start << "Cooperativity: " << _cooperativity << std::endl;
	
}


