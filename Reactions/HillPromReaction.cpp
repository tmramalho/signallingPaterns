/*
 *  HillPromReaction.cpp
 *	Signalling Patterns
 *
 *  Created by Michael Celentano 2 on 8/5/13.
 *
 */

#include "HillPromReaction.h"

/* Constructor: HillPromReaction(i_promoter,i_promoted,
 *								kinetic,K,cooperativity)
 * -------------------------------------------------------------------------- 
 * Constructs HillPromReaction, reading in data passed to it.
 */
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

HillPromReaction::HillPromReaction( std::ifstream& file ) {

	for (int i = 0 ; i < NUM_LINE ; i++) {
		switch (i) {
			case I_PROMOTER_LINE:
				file.ignore(256,':');
				file >> _i_promoter;
				break;
			case I_PROMOTED_LINE:
				file.ignore(256,':');
				file >> _i_promoted;
				break;
			case KINETIC_LINE:
				file.ignore(256,':');
				file >> _kinetic;
				break;
			case K_LINE:
				file.ignore(256,':');
				file >> _K;
				break;
			case COOPERATIVITY_LINE:
				file.ignore(256,':');
				file >> _cooperativity;
				break;
			default:
				break;
		}
	}
	
}

Reaction* HillPromReaction::copy() {
	return new HillPromReaction(_i_promoter,_i_promoted,_kinetic,_K,_cooperativity);
}

/* Public method: get_i_part(part_num)
 * -------------------------------------------------------------------------- 
 * Returns the index of the manager of the part_num participant in this
 * reaction.
 *
 *		Participant 0 is the promoted protein.
 *
 * Returns NEXIST for part_num != 0.
 */
int HillPromReaction::get_i_part( int part_num ) const {
	switch (part_num) {
		case 0:
			return _i_promoted;
			break;
		default:
			return NEXIST;
			break;
	}
}

int HillPromReaction::get_i_dependent_molecule() const {
	return NEXIST;
}

/* Public method: react(curr_tissue,dx_dt,i_curr_cell)
 * -------------------------------------------------------------------------- 
 * Adds to dx_dt vector (passed by reference), the deterministic rates of
 * change due to this PromBindingReaction, occuring in i_curr_cell.
 *		
 *		A promotes B
 *
 *		(d/dt)[B] = kinetic * ( [A]^n ) / ( K + [A]^n )
 *
 */
void HillPromReaction::react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell ) {

	double r = 
	_kinetic
	* pow( curr_tissue.at(i_curr_cell , _i_promoter ) , _cooperativity )
	/
	( _K + pow( curr_tissue.at(i_curr_cell , _i_promoter ) , _cooperativity ));
	 
	dx_dt.at(i_curr_cell,_i_promoted) += r;
}

/* Public method: react(curr_tissue,dx_dt,i_curr_cell
 *						dist,generator,q)
 * -------------------------------------------------------------------------- 
 * Adds to dx_dt vector (passed by reference), the stochastic rates of
 * change due to this PromBindingReaction, occuring in i_curr_cell.
 *		
 *		A promotes B
 *
 *		(d/dt)[B] = kinetic * ( [A]^n ) / ( K + [A]^n )
 *
 * Stochastic term is
 *
 *		det_flow * random_var * sqrt(q/dt)
 *
 * where q measures the magnitued of the noise.
 *
 * We divde by sqrt(dt) so that dx_dt * dt is the correct term for the 
 * the numerical stochastic integration of time step dt (see Documentation,
 * section INTEGRATION).
 */
void HillPromReaction::react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell ,
						 boost::random::normal_distribution<>& dist , 
						 boost::random::mt19937& generator , double q ) {
	
	/* dx = dt * f(x) + dt * (g(x) * rand * sqrt(q/dt)) */
	/* This implementation has g(x) = 1 */
	
	double det_flow = 
	_kinetic
	* pow( curr_tissue.at(i_curr_cell , _i_promoter ) , _cooperativity )
	/
	( _K + pow( curr_tissue.at(i_curr_cell , _i_promoter ) , _cooperativity ));	
	
	double rand = dist(generator);
	double stoc_flow = det_flow * rand * sqrt(q/(_sc_ref->_dt));
	
	double flow = det_flow + stoc_flow;
	
	if ( std::isnan(flow) ) {
		//std::cout << "Found NaN" << std::endl;
	}
	
	dx_dt.at(i_curr_cell,_i_promoted) += flow;
	
}

/* Public Method: update_indices(first_index,num_insertion)
 * -------------------------------------------------------------------------- 
 * Updates indices i_reac_zero, i_reac_one, i_product, given an insertion
 * of num_insertion molecules beginning at first_index into our genome.
 *
 * FOR NOW WE ONLY CONSIDER INSERTIONS, IE NUMINSERTIONS >= 0. 
 */
void HillPromReaction::update_mol_indices( int first_index , int num_insertion ) {
	/* Note: Because NEXIST = -1, non-existant participants will never be 
	 * updated by the insertion procedure as desired.
	 */		
	if (_i_promoter >= first_index) _i_promoter += num_insertion;
	if (_i_promoted >= first_index) _i_promoted += num_insertion;
}

/* Public Method: mutate(generator)
 * -------------------------------------------------------------------------- 
 * Mutates either the kinetic constant, K, or cooperativity, chosen
 * randomly. Mutates them through multiplication by random real chosen 
 * uniformly between 0 and 2.
 */
void HillPromReaction::mutate ( boost::random::mt19937& generator ) {

	boost::random::uniform_int_distribution<> which_kinetic(0,1);
	boost::random::uniform_real_distribution<> mutation_factor(0.0,2.0);
	
	int kinetic_num = which_kinetic(generator);
	double mut_factor = mutation_factor(generator);
	
	switch ( kinetic_num ) {
		case 0:
			_K *= mut_factor;
			break;
		case 1:
			_kinetic *= mut_factor;
			break;
		default:
			break;
	}
	
}

/* Public Method: print_info(line_start)
 * -------------------------------------------------------------------------- 
 * Prints the data describing the reaction, with line_start beginning each 
 * line.
 *
 * Format:
 *
 *		<line_start> Reaction Type: Hill Promotion Reaction
 *		<line_start> Index of Promoting Protein: <i>
 *		<line_start> Index of Promoted Protein: <i>
 *		<line_start> Kinetic Rate: <r>
 *		<line_start> K: <r>
 *		<line_start> Cooperativity: <r>
 *		
 */
void HillPromReaction::print_info ( std::string line_start ) {
	
	std::cout << line_start << "Reaction Type: Hill Promotion Reaction\n";
	
	for (int i = 0 ; i < NUM_LINE ; i++) {
		switch (i) {
			case I_PROMOTER_LINE:
				std::cout << line_start << "Index of Promoting Protein: " << _i_promoter << "\n";
				break;
			case I_PROMOTED_LINE:
				std::cout << line_start << "Index of Promoted Protein: " << _i_promoted <<  "\n";
				break;
			case KINETIC_LINE:
				std::cout << line_start << "Kinetic Rate: " << _kinetic <<  "\n";
				break;
			case K_LINE:
				std::cout << line_start << "K: " << _K << "\n";
				break;
			case COOPERATIVITY_LINE:
				std::cout << line_start << "Cooperativity: " << _cooperativity << "\n";
				break;
			default:
				break;
		}
	}
	
	
}

/* Public Method: to_file(file,line_start)
 * -------------------------------------------------------------------------- 	
 */
void HillPromReaction::to_file ( std::ofstream& file , std::string line_start ) {
	
	file << line_start << "Reaction Type: Hill Promotion Reaction\n";
	
	for (int i = 0 ; i < NUM_LINE ; i++) {
		switch (i) {
			case I_PROMOTER_LINE:
				file << line_start << "Index of Promoting Protein: " << _i_promoter << "\n";
				break;
			case I_PROMOTED_LINE:
				file << line_start << "Index of Promoted Protein: " << _i_promoted <<  "\n";
				break;
			case KINETIC_LINE:
				file << line_start << "Kinetic Rate: " << _kinetic <<  "\n";
				break;
			case K_LINE:
				file << line_start << "K: " << _K << "\n";
				break;
			case COOPERATIVITY_LINE:
				file << line_start << "Cooperativity: " << _cooperativity << "\n";
				break;
			default:
				break;
		}
	}
}

