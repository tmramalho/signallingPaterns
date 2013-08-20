/*
 *  PromReaction.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/22/13.
 *
 */

#include "PromReaction.h"

/* Constructor: PromReaction(i_gene,i_prot,kinetic)
 * -------------------------------------------------------------------------- 
 * Constructs PromReaction, reading in data passed to it.
 */
PromReaction::PromReaction( int i_gene , int i_prot , double kinetic ) {
	
	std::cout << "PromReaction Specific Constructor\n";
	
	_type = PROMOTION;
	
	_num_part = 1;
	
	_i_gene = i_gene;
	_i_prot = i_prot;
	_kinetic = kinetic;
	
}

PromReaction::~PromReaction() {
	std::cout << "PromReaction Destructor\n";
}

Reaction* PromReaction::copy() {
	return new PromReaction(_i_gene,_i_prot,_kinetic);
}

/* Public method: get_i_part(part_num)
 * -------------------------------------------------------------------------- 
 * Returns the index of the manager of the part_num participant in this
 * reaction.
 *
 *		Participant 0 is the protein being produced.
 *
 * Returns NEXIST for part_num != 0.
 */

int PromReaction::get_i_part( int part_num ) {
	switch (part_num) {
		case 0:
			return _i_prot;
			break;
		default:
			return NEXIST;
			break;
	}
}

/* Public method: react(curr_tissue,dx_dt,i_curr_cell)
 * -------------------------------------------------------------------------- 
 * Adds to dx_dt vector (passed by reference), the deterministic rates of
 * change due to this PromReaction, occuring in i_curr_cell.
 *		
 *		gene --> protein
 *		
 *		reaction: protein produced at rate kinetic*[gene]
 *
 */
void PromReaction::react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell ) {
	
	double r = _kinetic * curr_tissue.at( i_curr_cell , _i_gene );
	dx_dt.at( i_curr_cell , _i_prot ) += r;
	
}

/* Public method: react(curr_tissue,dx_dt,i_curr_cell
 *						dist,generator,q)
 * -------------------------------------------------------------------------- 
 * Adds to dx_dt vector (passed by reference), the stochastic rates of
 * change due to this PromReaction, occuring in i_curr_cell.
 *		
 *		gene --> protein
 *		
 *		reaction: protein produced at rate kinetic*[gene]
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
void PromReaction::react( dmat& curr_tissue , dmat& dx_dt , int i_curr_cell ,
						boost::random::normal_distribution<>& dist , 
						boost::random::mt19937& generator , double q ) {
	
	/* dx = dt * f(x) + dt * (g(x) * rand * sqrt(q/dt)) */
	/* This implementation has g(x) = 1 */
	
	double det_flow = _kinetic * curr_tissue.at( i_curr_cell , _i_gene ); 
	
	double rand = dist(generator);
	double stoc_flow = rand * sqrt(q/(_sc_ref->_dt));
	
	double flow = det_flow + stoc_flow;
	
	if ( std::isnan(flow) ) {
		//std::cout << "Found NaN" << std::endl;
	}
	
	dx_dt.at( i_curr_cell , _i_prot ) += flow;
	
}

/* Public Method: update_indices(first_index,num_insertion)
 * -------------------------------------------------------------------------- 
 * Updates indices i_reac_zero, i_reac_one, i_product, given an insertion
 * of num_insertion molecules beginning at first_index into our genome.
 *
 * FOR NOW WE ONLY CONSIDER INSERTIONS, IE NUMINSERTIONS >= 0. 
 */

void PromReaction::update_indices( int first_index , int num_insertion ) {
	/* Note: Because NEXIST = -1, non-existant participants will never be 
	 * updated by the insertion procedure as desired.
	 */	
	if (_i_gene >= first_index) _i_gene += num_insertion;
	if (_i_prot >= first_index) _i_prot += num_insertion;
}

/* Public Method: mutate(generator)
 * -------------------------------------------------------------------------- 
 * Mutates promotion rate through multiplication by random real chosen 
 * uniformly between 0 and 2.
 */
void PromReaction::mutate ( boost::random::mt19937& generator ) {
	
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
 *		<line_start> Reaction Type: Promotion Reaction
 *		<line_start> Index of Gene: <i>
 *		<line_start> Index of Produced Protein: <i>
 *		<line_start> Rate of Production: <r>
 *		
 */
void PromReaction::print_info ( std::string line_start ) {
	
	std::cout << line_start << "Reaction Type: Promotion Reaction" << std::endl;
	std::cout << line_start << "Index of Gene: " << _i_gene << std::endl;
	std::cout << line_start << "Index of Produced Protein: " << _i_prot << std::endl;
	std::cout << line_start << "Rate of Production: " << _kinetic << std::endl;
	
}

/* Public Method: to_file(file,line_start)
 * -------------------------------------------------------------------------- 
 */
void PromReaction::to_file ( std::ofstream& file , std::string line_start ) {

	file << line_start << "Reaction Type: Promotion Reaction\n";
	file << line_start << "Index of Gene: " << _i_gene << "\n";
	file << line_start << "Index of Produced Protein: " << _i_prot << "\n";
	file << line_start << "Rate of Production: " << _kinetic << "\n";
	
}

