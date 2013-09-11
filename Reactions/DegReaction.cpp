/*
 *  DegReaction.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/22/13.
 *
 */

#include "DegReaction.h"

/* Constructor: DegReaction(i_reac,kinetic)
 * -------------------------------------------------------------------------- 
 * Constructs DegReaction, reading in data passed to it.
 */
DegReaction::DegReaction( int i_reac , double kinetic ) {
	
	//std::cout << "DegReaction Specific Constructor\n";
	
	_type = DEGRADATION;
	
	_num_part = 1;
	
	_i_reac = i_reac;
	_kinetic = kinetic;
	
}

DegReaction::DegReaction( std::ifstream& file ) {

	for (int i = 0; i < NUM_LINE; i++) {
		switch (i) {
			case I_REAC_LINE:
				file.ignore(256,':');
				file >> _i_reac;
				break;
			case KINETIC_LINE:
				file.ignore(256,':');
				file >> _kinetic;
				break;
			default:
				break;
		}
	}

}

DegReaction::~DegReaction() {
	//std::cout << "DegReaction Destructor\n";
}

Reaction* DegReaction::copy() {
	return new DegReaction(_i_reac,_kinetic);
}

/* Public method: get_i_part(part_num)
 * -------------------------------------------------------------------------- 
 * Returns the index of the manager of the part_num participant in this
 * reaction.
 *
 *		Participant 0 is the protein being degraded.
 *
 * Returns NEXIST for part_num != 0.
 */
int DegReaction::get_i_part( int part_num ) const {
	switch (part_num) {
		case 0:
			return _i_reac;
			break;
		default:
			return NEXIST;
			break;
	}
}

int DegReaction::get_i_dependent_molecule() const {
	return NEXIST;
}

/* Public method: react(curr_tissue,dx_dt,i_curr_cell)
 * -------------------------------------------------------------------------- 
 * Adds to dx_dt vector (passed by reference), the deterministic rates of
 * change due to this DegReaction, occuring in i_curr_cell.
 *		
 *		reactant --> NOTHING
 *		
 *		forwardReaction: reactant degrades at rate -kinetic * [reactant]
 *
 */
void DegReaction::react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell ) {
	
	double r = _kinetic * curr_tissue.at(i_curr_cell,_i_reac);
	
	dx_dt.at(i_curr_cell,_i_reac) -= r;
	
}

/* Public method: react(curr_tissue,dx_dt,i_curr_cell,
 *						dist,generator,q)
 * -------------------------------------------------------------------------- 
 * Adds to dx_dt vector (passed by reference), the stochastic rates of
 * change due to this DegReaction, occuring in i_curr_cell.
 *		
 *		reactant --> NOTHING
 *		
 *		forwardReaction: reactant degrades at rate -kinetic * [reactant]
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
void DegReaction::react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell ,
						 boost::random::normal_distribution<>& dist , 
						 boost::random::mt19937& generator , double q ) {
	
	/* dx = dt * f(x) + dt * (g(x) * rand * sqrt(q/dt)) */
	/* This implementation has g(x) = det_flow */
	
	double det_flow = _kinetic * curr_tissue.at(i_curr_cell,_i_reac); 
	
	double rand = dist(generator);
	double stoc_flow = curr_tissue.at(i_curr_cell,_i_reac) * rand * sqrt(q/(_sc_ref->_dt));
	
	double flow = det_flow + stoc_flow;
	
	if ( std::isnan(flow) ) {
		//std::cout << "Found NaN" << std::endl;
	}
	
	dx_dt.at(i_curr_cell,_i_reac) -= flow;
	
}

/* Public Method: update_indices(first_index,num_insertion)
 * -------------------------------------------------------------------------- 
 * Updates indices i_reac_zero, i_reac_one, i_product, given an insertion
 * of num_insertion molecules beginning at first_index into our genome.
 *
 * FOR NOW WE ONLY CONSIDER INSERTIONS, IE NUMINSERTIONS >= 0. 
 */

void DegReaction::update_mol_indices( int first_index , int num_insertion ) {
	_i_reac = update_index( first_index , num_insertion , _i_reac );
}

/* Public Method: mutate(generator)
 * -------------------------------------------------------------------------- 
 * Mutates degredation kinetic constant through multiplication by random real 
 * chosen uniformly between 0 and 2.
 */
void DegReaction::mutate ( boost::random::mt19937& generator ) {
	
	boost::random::uniform_real_distribution<> mutation_factor(0.0,2.0);
	
	double mut_factor = mutation_factor(generator);
	
	_kinetic *= mut_factor;
	
}

/* Public Method: print_info(line_start)
 * -------------------------------------------------------------------------- 
 * Prints the data describing the reaction, with line_start beginning each 
 * line.
 *
 * Format:
 *
 *		<line_start> Reaction Type: Degredation Reaction
 *		<line_start> Index of Degrading Protein (in Complex): <i>
 *		<line_start> Degredation Rate: <r>
 *		
 */
void DegReaction::print_info ( std::string line_start ) {
	
	std::cout << line_start << "Reaction Type: Degradation Reaction\n";
	
	for (int i = 0; i < NUM_LINE; i++) {
		switch (i) {
			case I_REAC_LINE:
				std::cout << line_start << "Index of Degrading Protein: " << _i_reac << "\n";
				break;
			case KINETIC_LINE:
				std::cout << line_start << "Degredation Rate: " << _kinetic << "\n";
				break;
			default:
				break;
		}
	}
	
}

/* Public Method: to_file(file,line_start)
 * -------------------------------------------------------------------------- 
 */
void DegReaction::to_file ( std::ofstream& file , std::string line_start ) {
	
	file << line_start << "Reaction Type: Degradation Reaction\n";
	
	for (int i = 0; i < NUM_LINE; i++) {
		switch (i) {
			case I_REAC_LINE:
				file << line_start << "Index of Degrading Protein: " << _i_reac << "\n";
				break;
			case KINETIC_LINE:
				file << line_start << "Degredation Rate: " << _kinetic << "\n";
				break;
			default:
				break;
		}
	}
}

