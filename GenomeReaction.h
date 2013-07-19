/*
 *  GenomeReaction.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/18/13.
 *
 */



# include "ReactionType.h"
# include "MoleculeType.h"

#ifndef GENOMEREACTION_H
#define GENOMEREACTION_H

class GenomeReaction {
	
public:
	GenomeReaction();
	~GenomeReaction();
	
	int getNumPart();
	int getIPart(int partNum);
	ReactionType getReacType();
	MoleculeType getMolType(int partNum);
	
	double getKinetic(int num);
	
private:
	
	static const int NEXIST = -1;
	
	// ODEReaction will react differently depending on the reaction type.
	ReactionType type;
	
	// There are between one and three participants in a reaction.
	int numPart;
	
	// We will assign NEXIST to the location of any non-existant participant.
	int iPartZero;
	MoleculeType typePartZero;
	
	int iPartOne;
	MoleculeType typePartOne;

	int iPartTwo;
	MoleculeType typePartTwo;
	
	double kineticZero;
	double kineticOne;
	
};

#endif