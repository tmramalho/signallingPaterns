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
Molecule::Molecule() {}

/* Destructor: ~Molecule()
 * -------------------------------------------------------------------------- 
 */
Molecule::~Molecule() {}

/* Public Method: getI()
 * -------------------------------------------------------------------------- 
 * Returns the index of this molecule in the genome that contains it. Because
 * the genome contains separate vectors for genes and proteins, the index i
 * refers only to its index in the vector that contains it. To know in which
 * vector to look, we also need its type.
 */
unsigned int Molecule::getISelf() {
	return iSelf;
}

/* Public Method: getInitConc()
 * -------------------------------------------------------------------------- 
 * Returns the initial concentration of the gene in the simulation.
 */

double Molecule::getInitConc() {
	return initConc;
}

/* Public Method: addReaction(reacRef)
 * -------------------------------------------------------------------------- 
 * Adds to the list the molecule keeps of the reactions it participates in a
 * reaction. 
 *
 * CAUTION: The molecule does not own the reaction to which it points. Rather,
 * the genome does. The molecule, even in its destructor, should never delete
 * these reactions.
 */
void Molecule::addReaction(Reaction* const reacRef){
	reactions.push_back(reacRef);
}


