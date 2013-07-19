/*
 *  Molecule.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/17/13.
 *
 */

#include "Molecule.h"

Molecule::Molecule() {}

Molecule::~Molecule() {}

/* Public Method: getLocation()
 * -------------------------------------------------------------------------- 
 * Returns the location of this molecule in either the genes or proteins 
 * genome of this molecule.
 */

unsigned int Molecule::getI() {
	return i;
}

void Molecule::addReaction(GenomeReaction* const reactionRef){
	reactions.push_back(reactionRef);
}


