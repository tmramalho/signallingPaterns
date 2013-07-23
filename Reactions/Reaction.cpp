/*
 *  Reaction.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/18/13.
 *
 */

#include "Reaction.h"


Reaction::Reaction() {}

Reaction::~Reaction() {}

/* Public Method: getNumPart()
 * --------------------------------------------------------------------------
 * Returns the number of participants in this reaction. 
 */
int Reaction::getNumPart() {
	return numPart;
}
