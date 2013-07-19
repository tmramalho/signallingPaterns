/*
 *  GenomeReaction.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/18/13.
 *
 */

# include "GenomeReaction.h"

GenomeReaction::GenomeReaction() {
	type = COMBINATION;
	
	numPart = 3;
	
	iPartZero = 0;
	typePartZero = PROTEIN;
	
	iPartOne = 1;
	typePartOne = PROTEIN;
	
	iPartTwo = 2;
	typePartTwo = PROTEIN;
	
	kineticZero = 0.7;
	kineticOne = 0.3;
}


GenomeReaction::~GenomeReaction() {}

/* Public Method: getNumPart()
 * -------------------------------------------------------------------------- 
 * Returns the number of participating molecules in this reaction.
 */
int GenomeReaction::getNumPart() {
	return numPart;
}

/* Public Method: getMolLoc(partNum)
 * -------------------------------------------------------------------------- 
 * Returns the location in the vector of molecules in the genome (either the
 * vector of proteins or of genes) of the participant in the reaction indexed 
 * by numParticipant.
 */
int GenomeReaction::getIPart(int partNum) {
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

/* Public Method: getReacType()
 * -------------------------------------------------------------------------- 
 * Returns a ReactionType describing what type of reaction this is.
 */
ReactionType GenomeReaction::getReacType() {
	return type;
}

/* Public Method: getMolType(partNum)
 * -------------------------------------------------------------------------- 
 * Returns a MoleculeType describing whether this participant is a gene or
 * a protein.
 */
MoleculeType GenomeReaction::getMolType(int partNum) {
	switch (partNum) {
		case 0:
			return typePartZero;
			break;
		case 1:
			return typePartOne;
			break;
		case 2:
			return typePartTwo;
			break;
	}
}

/* Public Method: returnKinetic(num)
 * -------------------------------------------------------------------------- 
 * Returns the value of the num kinetic constant for this reaction.
 */

double GenomeReaction::getKinetic(int num) {
	switch (num) {
		case 0:
			return kineticZero;
			break;
		case 1:
			return kineticOne;
			break;
	}
}


