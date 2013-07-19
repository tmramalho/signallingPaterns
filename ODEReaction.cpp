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
 *		
 *		moleculeZero + moleculeOne <--> moleculeTwo
 *		forwardRate = 0.7 * [moleculeZero] * [moleculeOne]
 *		backwardRate = 0.3 * [moleculeTwo]
 *
 */
ODEReaction::ODEReaction() {
	
	type = COMBINATION;
	
	numPart = 3;
	
	// The first cell has index zero, in which all participants are in.
	cellPartZero = 0;
	cellPartOne = 0;
	cellPartTwo = 0;
	
	iPartZero = 0;
	iPartOne = 1;
	iPartTwo = 2;
	
	dxdtPartZero = 0.0;
	dxdtPartOne = 0.0;
	dxdtPartTwo = 0.0;
	
	kineticZero = 0.7;
	kineticOne = 0.3;
	
}

/* Constructor: ODEReaction(---full argument list---)
 * --------------------------------------------------------------------------
 * Brute force constructor that simply reads in all the data the reaction
 * needs.
 */
ODEReaction::ODEReaction( ReactionType type , 
			int numParticipants , 
			int cellPartZero , int cellPartOne , int cellPartTwo ,
			int iPartZero , int iPartOne , int iPartTwo ,
			int dxdtPartZero , int dxdtPartOne , int dxdtPartTwo , 
			double kineticZero , double kineticOne ) {
	
	this->type = type;
	
	this->cellPartZero = cellPartZero;
	this->cellPartOne = cellPartOne;
	this->cellPartTwo = cellPartTwo;
	
	this->iPartZero = iPartZero;
	this->iPartOne = iPartOne;
	this->iPartTwo = iPartTwo;
	
	this->dxdtPartZero = dxdtPartZero;
	this->dxdtPartOne = dxdtPartOne;
	this->dxdtPartTwo = dxdtPartTwo;
	
	this->kineticZero = kineticZero;
	this->kineticOne = kineticOne;
	
}

/* Constructor: ODEReaction(genome,iReac,iCell)
 * --------------------------------------------------------------------------
 * This constructor creates the ODEReaction corresponding to the GenomeReaction
 * held in index iReac of genome and which operates in the cell with index
 * iCell.
 */
ODEReaction::ODEReaction ( Genome& genome , int iReac , int iCell ) {
	
	GenomeReaction* reacRef = genome.getReacRef(iReac);
	
	type = reacRef->getReacType();
	
	numPart = reacRef->getNumPart();
	
	iPartZero = genome.getIPartForODE(iReac,0);
	iPartOne = genome.getIPartForODE(iReac,1);
	iPartTwo = genome.getIPartForODE(iReac,2);
	
	cellPartZero = iCell;
	cellPartOne = iCell;
	cellPartTwo = iCell;
	
	kineticZero = reacRef->getKinetic(0);
	kineticOne = reacRef->getKinetic(1);
	
}

/* Destructor
 * --------------------------------------------------------------------------
 */
ODEReaction::~ODEReaction() {}

/* Public Method: getNumPart()
 * --------------------------------------------------------------------------
 * Returns the number of participants in this reaction. It should only be 
 * 1, 2, or 3. 
 */
int ODEReaction::getNumPart() {
	return numPart;
}

// We assume we never call this for a numParticipant above numParticipants
int ODEReaction::getICell(int partNum) {
	switch (partNum) {
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

int ODEReaction::getIPart(int partNum) {
	switch (partNum) {
		case 0:
			return iPartZero;
			break;
		case 1:
			return iPartOne;
			break;
		case 2:
			return iPartTwo;
			break;
	}
}

// If this is called before we call the react method, it may not return the 
// appropriate current change rate.
double ODEReaction::getDxDt(int partNum) {
	switch (partNum) {
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
			currTissue.at(cellPartZero)->at(iPartZero) *
			currTissue.at(cellPartOne)->at(iPartOne);
			
			double backwardFlow = kineticOne *
			currTissue.at(cellPartTwo)->at(iPartTwo);

			/* Update our dxdt values */
			dxdtPartZero = dxdtPartOne = backwardFlow - forwardFlow;
			dxdtPartTwo = forwardFlow - backwardFlow;
			
			break;
	}
}
















