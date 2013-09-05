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

Gene::Gene( const Gene& other ) : Molecule(other) {
	_i_product = other._i_product;
	_i_bound_promoter = other._i_bound_promoter;
	_i_root = other._i_root;
}

Gene::Gene( std::ifstream& file ) {
	
	for (int i = 0; i < NUM_LINE ; i++) {
		switch (i) {
			case I_SELF_LINE:
				file.ignore(256,':');
				file >> _i_self;
				break;
			case I_PRODUCT_LINE:
				file.ignore(256,':');
				file >> _i_product;
				break;
			case I_BOUND_PROMOTER_LINE:
				file.ignore(256,':');
				file >> _i_bound_promoter;
				break;
			case I_ROOT_LINE:	
				file.ignore(256,':');
				file >> _i_root;
				break;
			case INIT_CONC_LINE:
				file.ignore(256,':');
				file >> _init_conc;
				break;
			default:
				break;
		}
	}
	
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

void Gene::update_mol_indices( int first_index , int num_insertion ) {
	
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
	for (int i = 0; i < NUM_LINE ; i++) {
		switch (i) {
			case I_SELF_LINE:
				std::cout << line_start << "Index of Self: " << _i_self << "\n";
				break;
			case I_PRODUCT_LINE:
				std::cout << line_start << "Index of Protein This Produces: " << _i_product << "\n";
				break;
			case I_BOUND_PROMOTER_LINE:
				std::cout << line_start << "Index of Protein Bound to This Gene: " << _i_bound_promoter << "\n";
				break;
			case I_ROOT_LINE:	
				std::cout << line_start << "Index of Gene This Gene Derived From: " << _i_root << "\n";
				break;
			case INIT_CONC_LINE:
				std::cout << line_start << "Initial Concentration: " << _init_conc << "\n";
				break;
			default:
				break;
		}
	}
}

/* Public Method: to_file(file)
 * -------------------------------------------------------------------------- 
 */

void Gene::to_file( std::ofstream& file , std::string line_start ) {
	
	for (int i = 0; i < NUM_LINE ; i++) {
		switch (i) {
			case I_SELF_LINE:
				file << line_start << "Index of Self: " << _i_self << "\n";
				break;
			case I_PRODUCT_LINE:
				file << line_start << "Index of Protein This Produces: " << _i_product << "\n";
				break;
			case I_BOUND_PROMOTER_LINE:
				file << line_start << "Index of Protein Bound to This Gene: " << _i_bound_promoter << "\n";
				break;
			case I_ROOT_LINE:	
				file << line_start << "Index of Gene This Gene Derived From: " << _i_root << "\n";
				break;
			case INIT_CONC_LINE:
				file << line_start << "Initial Concentration: " << _init_conc << "\n";
				break;
			default:
				break;
		}
	}
	
}




