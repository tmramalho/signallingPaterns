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
	
	/* For reactions with less than three participants, all indexes and 
	 * locations for out-of-scope participants will be set to NEXIST.
	 */
	static const int NEXIST = -1;
	
	/* ODEReaction will calculate rates differently depending on its type */
	ReactionType type;
	
	/* Number of Participants: must be betwee 1 and 3 */
	int numPart;
	
	/* The type and index for each participant.
	 * In the genome, genes are held in a different vector than proteins. The
	 * index is its index in whichever vector it is in (Thus, presumably, we
	 * could have a situation where iPartZero = iPartOne if these participants
	 * are of different types.) Which vector to look in is determined by the
	 * MoleculeType.
	 */
	int iPartZero;
	MoleculeType typePartZero;
	
	int iPartOne;
	MoleculeType typePartOne;

	int iPartTwo;
	MoleculeType typePartTwo;
	
	/* Kinetic constants for the reaction. Depending on the type of reaction,
	 * these may have different interpretations.
	 */
	double kineticZero;
	double kineticOne;
	
};

#endif