/*
 *  Molecule.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/17/13.
 *
 */

#include "Molecule.h"

/* Constructor: Molecule()
 * -------------------------------------------------------------------------- 
 */
Molecule::Molecule() {
	_sc_ref = SettingsCont::getInstance();
}

/* Public Method: get_i_self()
 * -------------------------------------------------------------------------- 
 * Returns the index of this molecule in the genome that contains it. Because
 * the genome contains separate vectors for genes and proteins, the index i
 * refers only to its index in the vector that contains it. To know in which
 * vector to look, we also need its type.
 */
unsigned int Molecule::get_i_self() {
	return _i_self;
}

/* Public Method: get_init_conc()
 * -------------------------------------------------------------------------- 
 * Returns the initial concentration of the gene in the simulation.
 */

double Molecule::get_init_conc() {
	return _init_conc;
}

/* Public Method: add_reaction(reac_ref)
 * -------------------------------------------------------------------------- 
 * Adds to the list the molecule keeps of the reactions it participates in a
 * reaction. 
 *
 * CAUTION: The molecule does not own the reaction to which it points. 
 * Rather, the manager does. The molecule, even in its destructor, should 
 * never delete these reactions.
 */
void Molecule::add_reaction(Reaction* reac_ref){
	_reactions.push_back(reac_ref);
}

/* Public Method: get_reaction(i_reac)
 * -------------------------------------------------------------------------- 
 * Returns a reference to the reaction stored at index i_reac of the 
 * molecule.
 *
 * We assume i_reac is a valid index. No error checking is performed.
 */
Reaction* Molecule::get_reaction( int i_reac ) {
	return _reactions.at( i_reac );
}

/* Public Method: get_num_reac()
 * -------------------------------------------------------------------------- 
 * Returns the number of reactions the molecule participates in. This is 
 * useful for avoiding passing an invalid index to get_reaction.
 */
int Molecule::get_num_reac() {
	return _reactions.size();
}






