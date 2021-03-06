/*
 *  HillDegReaction.cpp
 *  Signalling Patterns
 *
 *  Created by Michael Celentano 2 on 8/5/13.
 *
 */

#include "HillRepReaction.h"

/* Constructor: HillPromReaction(i_repressor,i_repressed,
 *								kinetic,K,cooperativity)
 * -------------------------------------------------------------------------- 
 * Constructs HillPromReaction, reading in data passed to it.
 */
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

HillRepReaction::HillRepReaction( std::ifstream& file ) {

	for (int i = 0 ; i < NUM_LINE ; i++) {
		switch (i) {
			case I_REPRESSOR_LINE:
				file.ignore(256,':');
				file >> _i_repressor;
				break;
			case I_REPRESSED_LINE:
				file.ignore(256,':');
				file >> _i_repressed;
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

Reaction* HillRepReaction::copy() {
	return new HillRepReaction(_i_repressor,_i_repressed,_kinetic,_K,_cooperativity);
}

/* Public method: get_i_part(part_num)
 * -------------------------------------------------------------------------- 
 * Returns the index of the manager of the part_num participant in this
 * reaction.
 *
 *		Participant 0 is the repressed protein.
 *
 * Returns NEXIST for part_num != 0.
 */
int HillRepReaction::get_i_part( int part_num ) const {
	switch (part_num) {
		case 0:
			return _i_repressed;
			break;
		default:
			return NEXIST;
			break;
	}
}

int HillRepReaction::get_i_dependent_molecule() const {
	return NEXIST;
}

/* Public method: react(curr_tissue,dx_dt,i_curr_cell)
 * -------------------------------------------------------------------------- 
 * Adds to dx_dt vector (passed by reference), the deterministic rates of
 * change due to this PromBindingReaction, occuring in i_curr_cell.
 *		
 *		A repressed B
 *
 *		(d/dt)[B] = kinetic * K / ( K + [A]^n )
 *
 */
void HillRepReaction::react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell ) {
	
	double r = 
	_kinetic
	* _K
	/
	( _K + pow( curr_tissue.at(i_curr_cell , _i_repressor ) , _cooperativity ));
	
	dx_dt.at(i_curr_cell,_i_repressed) += r;
}

/* Public method: react(curr_tissue,dx_dt,i_curr_cell
 *						dist,generator,q)
 * -------------------------------------------------------------------------- 
 * Adds to dx_dt vector (passed by reference), the stochastic rates of
 * change due to this PromBindingReaction, occuring in i_curr_cell.
 *		
 *		A repressed B
 *
 *		(d/dt)[B] = kinetic * K / ( K + [A]^n )
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
void HillRepReaction::react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell ,
						 boost::random::normal_distribution<>& dist , 
						 boost::random::mt19937& generator , double q ) {
	
	/* dx = dt * f(x) + dt * (g(x) * rand * sqrt(q/dt)) */
	/* This implementation has g(x) = 1 */
	
	double det_flow = 
	_kinetic
	* _K
	/
	( _K + pow( curr_tissue.at(i_curr_cell , _i_repressor ) , _cooperativity ));
	
	double rand = dist(generator);
	double stoc_flow = curr_tissue.at(i_curr_cell,_i_repressed) * rand * sqrt(q/(_sc_ref->_dt));
	
	double flow = det_flow + stoc_flow;
	
	if ( std::isnan(flow) ) {
		//std::cout << "Found NaN" << std::endl;
	}
	
	dx_dt.at(i_curr_cell,_i_repressed) += flow;
	
}

/* Public Method: update_indices(first_index,num_insertion)
 * -------------------------------------------------------------------------- 
 * Updates indices i_reac_zero, i_reac_one, i_product, given an insertion
 * of num_insertion molecules beginning at first_index into our genome.
 *
 * FOR NOW WE ONLY CONSIDER INSERTIONS, IE NUMINSERTIONS >= 0. 
 */
void HillRepReaction::update_mol_indices( int first_index , int num_insertion ) {
	_i_repressor = update_index( first_index , num_insertion , _i_repressor );
	_i_repressed = update_index( first_index , num_insertion , _i_repressed );
}

/* Public Method: mutate(generator)
 * -------------------------------------------------------------------------- 
 * Mutates either the kinetic constant, K, or cooperativity, chosen
 * randomly. Mutates them through multiplication by random real chosen 
 * uniformly between 0 and 2.
 */
void HillRepReaction::mutate ( boost::random::mt19937& generator ) {
	
	boost::random::uniform_int_distribution<> which_kinetic(0,1);
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
 *		<line_start> Reaction Type: Hill Repression Reaction
 *		<line_start> Index of Repressor Protein: <i>
 *		<line_start> Index of Repressed Protein: <i>
 *		<line_start> Kinetic Rate: <r>
 *		<line_start> K: <r>
 *		<line_start> Cooperativity: <r>
 *		
 */
void HillRepReaction::print_info ( std::string line_start ) {
	
	std::cout << line_start << "ReactionType: Hill Repression Reaction\n";
	
	for (int i = 0 ; i < NUM_LINE ; i++) {
		switch (i) {
			case I_REPRESSOR_LINE:
				std::cout << line_start << "Index of Repressor Protein: " << _i_repressor << "\n";
				break;
			case I_REPRESSED_LINE:
				std::cout << line_start << "Index of Repressed Protein: " << _i_repressed << "\n";
				break;
			case KINETIC_LINE:
				std::cout << line_start << "Kinetic Rate: " << _kinetic << "\n";
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
void HillRepReaction::to_file ( std::ofstream& file , std::string line_start ) {
	
	file << line_start << "ReactionType: Hill Repression Reaction\n";
	
	for (int i = 0 ; i < NUM_LINE ; i++) {
		switch (i) {
			case I_REPRESSOR_LINE:
				file << line_start << "Index of Repressor Protein: " << _i_repressor << "\n";
				break;
			case I_REPRESSED_LINE:
				file << line_start << "Index of Repressed Protein: " << _i_repressed << "\n";
				break;
			case KINETIC_LINE:
				file << line_start << "Kinetic Rate: " << _kinetic << "\n";
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



