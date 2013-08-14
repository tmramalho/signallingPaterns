/*
 *  Gene.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/17/13.
 *
 */

#include "Gene.h"

/* Constructor: Gene(i_self,i_product,i_bound_promoter,i_root,init_conc)
 * -------------------------------------------------------------------------- 
 * Constructs gene, reading in data passed to it.
 */
Gene::Gene( int i_self, int i_product , int i_bound_promoter , int i_root , double init_conc ) {
	
	
	_i_product = i_product;
	_i_bound_promoter = i_bound_promoter;
	_i_root = i_root;
	_i_self = i_self;
	_init_conc = init_conc;
}

/* Public Method: get_i_product()
 * -------------------------------------------------------------------------- 
 * Returns the index of this gene in the dvecs held in the cell.
 */
int Gene::get_i_product() {
	return _i_product;
}

/* Public Method: get_i_bound_promoter()
 * -------------------------------------------------------------------------- 
 * Returns the index in the dvecs of the protein that our bound to this 
 * gene as a promoter (or represser). If there is no such protein, returns
 * NEXIST.
 */
int Gene::get_i_bound_promoter() {
	return _i_bound_promoter;
}

/* Public Method: get_i_root()
 * -------------------------------------------------------------------------- 
 * Returns the index in the dvecs of the gene from which this gene derives,
 * if this gene is complexed with a protein. If it is not, it returns 
 * NEXIST.
 *
 * Example: For the gene-protein complex a:B, this would return the index of
 * a
 */
int Gene::get_i_root() {
	return _i_root;
}
 
/* Public Method: print_info(line_start)
 * -------------------------------------------------------------------------- 
 * Prints the data describing the gene, with line_start beginning each line.
 *
 * Format:
 *
 *		<line_start> Index of Self: <i>
 *		<line_start> Index of Protein Bound to This Gene: <i>
 *		<line_start> Index of Gene This Gene Derived From: <i>
 *		
 */

void Gene::print_info( std::string line_start ) {
	std::cout << line_start << "Index of Self: " << _i_self << std::endl;
	std::cout << line_start << "Index of Protein This Produces: " << _i_product << std::endl;
	if ( _i_bound_promoter != NEXIST ) {
		std::cout << line_start << "Index of Protein Bound to This Gene: " << _i_bound_promoter << std::endl;
	}
	else {
		std::cout << line_start << "Index of Protein Bound to This Gene: NO PROTEINS BOUND TO THIS GENE" << std::endl;
	}
	if ( _i_root != NEXIST ) {
		std::cout << line_start << "Index of Gene This Gene Derived From: " << _i_root << std::endl;
	}
	else {
		std::cout << line_start << "Index of Gene This Gene Derived From: GENE NOT DERIVED FROM SIMPLER GENE" << std::endl;
	}
}

/* Public Method: to_file(file)
 * -------------------------------------------------------------------------- 
 */

void Gene::to_file( std::ofstream& file , std::string line_start ) {
	file << line_start << "Index of Self: " << _i_self << "\n";
	file << line_start << "Index of Protein This Produces: " << _i_product << "\n";
	if ( _i_bound_promoter != NEXIST ) {
		file << line_start << "Index of Protein Bound to This Gene: " << _i_bound_promoter << "\n";
	}
	else {
		file << line_start << "Index of Protein Bound to This Gene: NO PROTEINS BOUND TO THIS GENE" << "\n";
	}
	if ( _i_root != NEXIST ) {
		file << line_start << "Index of Gene This Gene Derived From: " << _i_root << "\n";
	}
	else {
		file << line_start << "Index of Gene This Gene Derived From: GENE NOT DERIVED FROM SIMPLER GENE" << "\n";
	}
}

/* Private Method: update_indices(first_index,num_insertion)
 * -------------------------------------------------------------------------- 
 * Updates the index of i_self, i_product, i_bound_promoter, i_root, given an 
 * insertion of size num_insertion of molecules beginning at first_index.
 *
 * FOR NOW WE ONLY CONSIDER INSERTIONS, IE NUMINSERTIONS >= 0. 
 */

void Gene::update_indices( int first_index , int num_insertion ) {
	/* Because NEXIST = -1, it will never be updated by our insertion 
	 * procedure, as desired. */
	if ( _i_self >= first_index ) _i_self += num_insertion;
	if ( _i_product >= first_index ) _i_product += num_insertion;
	if ( _i_bound_promoter >= first_index ) _i_bound_promoter += num_insertion;
	if ( _i_root >= first_index ) _i_root += num_insertion;
}


