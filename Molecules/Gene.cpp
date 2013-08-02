/*
 *  Gene.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/17/13.
 *
 */

#include "Gene.h"

Gene::Gene( int i_self, int i_product , int i_bound_promoter , int i_root , double init_conc ) {
	_i_product = i_product;
	_i_bound_promoter = i_bound_promoter;
	_i_root = i_root;
	_i_self = i_self;
	_init_conc = init_conc;
}

/* Destructor: ~Gene()
 * -------------------------------------------------------------------------- 
 * No heap allocated memory owned by the gene class.
 */
Gene::~Gene() {}

/* Public Method: getIProduct()
 * -------------------------------------------------------------------------- 
 * Returns the index of this gene in the dvecs held in the cell.
 */
int Gene::get_i_product() {
	return _i_product;
}

/* Public Method: getIBoundPromoter()
 * -------------------------------------------------------------------------- 
 * Returns the index in the dvecs of the protein that our bound to this 
 * gene as a promoter (or represser). If there is no such protein, returns
 * NEXIST.
 */
int Gene::get_i_bound_promoter() {
	return _i_bound_promoter;
}

/* Public Method: getIRoot()
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

/* Private Method: updateIndices(first_index,num_insertion)
 * -------------------------------------------------------------------------- 
 * Updates the index of iSelf, iProduct, iBoundPromoter, iRoot, given an 
 * insertion of size numInsertions, beginning at firstIndex, into our dvecs 
 * containing molecule concentrations in our manager.
 *
 * See description for updateIndices private method of the manager class in
 * the Manager.cpp file for more precise description of insertion process.
 *
 * Note: Because NEXIST = -1, it will never be updated by insertion procedure,
 * as desired.
 *
 * FOR NOW WE ONLY CONSIDER INSERTIONS, IE NUMINSERTIONS >= 0. 
 */

void Gene::update_indices( int first_index , int num_insertion ) {
	if ( _i_self >= first_index ) _i_self += num_insertion;
	if ( _i_product >= first_index ) _i_product += num_insertion;
	if ( _i_bound_promoter >= first_index ) _i_bound_promoter += num_insertion;
	if ( _i_root >= first_index ) _i_root += num_insertion;
}


