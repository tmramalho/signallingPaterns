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
	
	/* The first cell has index zero, in which all participants are in. */
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
 * This class allocates nothing on the heap.
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

/* Public Method: getICell(partNum)
 * --------------------------------------------------------------------------
 * Returns the index of the cell of the partNum participant in this reaction.
 *
 * CAUTION: We assume we never call this for a partNum above numPart (the 
 * number of participants in the reaction). If we do, this will return
 * NEXIST. If partNum is not among 1, 2, or 3, nothing returns.
 */
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

/* Public Method: getIPart(partNum)
 * --------------------------------------------------------------------------
 * Returns the index in the dvec held by ODEManager of the partNum 
 * participant in this reaction.
 *
 * For example: If our cells contain four substances [A, B, C, and D] and
 * there is a COMBINATION reaction A + D <--> B, getIPart(2) == 3 because
 * D is the second participant in the COMBINATION REACTION and is in 
 * location 3 of the dvec held by the cell.
 *
 * CAUTION: We assume we never call this for a partNum above numPart (the 
 * number of participants in the reaction). If we do, this will return
 * NEXIST. If partNum is not among 1, 2, or 3, nothing returns.
 */
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

/* Public Method: getDxDt(partNum)
 * --------------------------------------------------------------------------
 * Returns the current value of dxdt for the concentration of the partNum 
 * participant in the reaction.
 *
 * CAUTION: No calculation is performed, so if the state of our cells has
 * changed since the last time we asked ODEReaction to react, this will not
 * return the current value of dxdt. It is best to always call react before
 * calling this method.
 */
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
 * Updates the dxdt values held by the ODEReaction based on the 
 * concentrations in the currTissue it receives. Nothing is changed in the
 * currTissue vector, nor in the dxdt vector of the ODEManager which contains
 * contains it.
 */
void ODEReaction::react( const std::vector< dvec* >& currTissue ) {
	
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
















