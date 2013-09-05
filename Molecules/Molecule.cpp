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
Molecule::Molecule() {_sc_ref = SettingsCont::getInstance();}

Molecule::Molecule( const Molecule& other ) {
	_sc_ref = SettingsCont::getInstance();
	_i_self = other._i_self;
	_init_conc = other._init_conc;
	_reactions = other._reactions;
}

/* Public Method: get_i_self()
 * -------------------------------------------------------------------------- 
 * Returns the index of this molecule in the genome that contains it. Because
 * the genome contains separate vectors for genes and proteins, the index i
 * refers only to its index in the vector that contains it. To know in which
 * vector to look, we also need its type.
 */
int Molecule::get_i_self() {return _i_self;}

/* Public Method: get_init_conc()
 * -------------------------------------------------------------------------- 
 * Returns the initial concentration of the gene in the simulation.
 */

double Molecule::get_init_conc() {return _init_conc;}

/* Public Method: add_reaction(reac_ref)
 * -------------------------------------------------------------------------- 
 * Adds to the list the molecule keeps of the reactions it participates in a
 * reaction. 
 *
 * CAUTION: The molecule does not own the reaction to which it points. 
 * Rather, the manager does. The molecule, even in its destructor, should 
 * never delete these reactions.
 */
void Molecule::add_reaction( int i_reac ) {_reactions.insert(i_reac);}

template <typename UpdateClass>
void Molecule::update_reac_indices( UpdateClass update_index ) {
	
	std::set<int> updates;
	
	for (std::set<int>::iterator it = _reactions.begin(); it != _reactions.end(); it++) {
		updates.insert(update_index(*it));
	}
	
	_reactions = updates;
	_reactions.erase(NEXIST);	
}

std::set<int>::iterator Molecule::reacs_begin() const {return _reactions.begin();}

std::set<int>::iterator Molecule::reacs_end() const {return _reactions.end();}
