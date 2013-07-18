/*
 *  ODEReaction.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/18/13.
 *
 */

#include "ODEReaction.h"


/* Constructor: ODEReaction()
 * --------------------------------------------------------------------------
 * This default constructor will create a COMBINATION reaction between the
 * three molecules in the first cell of a tissue.
 */
ODEReaction::ODEReaction() {
	
	type = COMBINATION;
	
	numParticipants = 3;
	
	// The first cell has index zero, in which all participants are in.
	cellPartZero = 0;
	cellPartOne = 0;
	cellPartTwo = 0;
	
	locPartZero = 0;
	locPartOne = 1;
	locPartTwo = 2;
	
	dxdtPartZero = 0.0;
	dxdtPartOne = 0.0;
	dxdtPartTwo = 0.0;
	
	kineticZero = 0.7;
	kineticOne = 0.3;
	
}

ODEReaction::~ODEReaction() {}

int ODEReaction::getNumPart() {
	return numParticipants;
}

// We assume we never call this for a numParticipant above numParticipants
int ODEReaction::getCellLoc(int numParticipant) {
	switch (numParticipant) {
		case 0:
			return cellPartZero;
			break;
		case 1:
			return cellPartOne;
			break;
		case 2:
			return cellPartTwo;
			break;
	}
}

int ODEReaction::getMolLoc(int numParticipant) {
	switch (numParticipant) {
		case 0:
			return locPartZero;
			break;
		case 1:
			return locPartOne;
			break;
		case 2:
			return locPartTwo;
			break;
	}
}

// If this is called before we call the react method, it may not return the 
// appropriate current change rate.
double ODEReaction::getDxDt(int numParticipant) {
	switch (numParticipant) {
		case 0:
			return dxdtPartZero;
			break;
		case 1:
			return dxdtPartOne;
			break;
		case 2:
			return dxdtPartTwo;
			break;
	}
}


/* Public Method: react(currTissue)
 * --------------------------------------------------------------------------
 */
void ODEReaction::react( std::vector< dvec* >& currTissue ) {
	
	// Might there be a faster implementation of this using inheritance.
	// For now let's do this to keep it simpler.
	switch (type) {
			
		case PROMOTION:
			break;
			
		case DEGRADATION:
			break;
			
		case COMBINATION:
			
			/* These reactions are of the form:
			 *
			 *		participantZero + participantOne <--> participantTwo
			 *		
			 *		forwardReaction: participantZero and participantOne combine 
			 *		to form participantTwo with rate kineticZero * 
			 *		[participantOne] * [particpantTwo].
			 *		
			 *		backwardReaction: participantTwo degrades into
			 *		participantZero and participantOne with rate kineticOne *
			 *		[participantTwo].
			 *
			 */
			
			/* Find forward and backward flows */
			double forwardFlow = kineticZero * 
			currTissue.at(cellPartZero)->at(locPartZero) *
			currTissue.at(cellPartOne)->at(locPartOne);
			
			double backwardFlow = kineticOne *
			currTissue.at(cellPartTwo)->at(locPartTwo);

			/* Update our dxdt values */
			dxdtPartZero = dxdtPartOne = backwardFlow - forwardFlow;
			dxdtPartTwo = forwardFlow - backwardFlow;
			
			break;
	}
}
















