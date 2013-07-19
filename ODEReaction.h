/*
 *  ODEReaction.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/18/13.
 *
 */

/* ODEReaction class
 * -------------------------------------------------------------------------- 
 * The ODEReaction contains information about where to find all of its
 * participants in the ODEManager, and a function to then calculate the rates
 * of change in concentration of all of its participants.
 *
 * The main purpose of the class is to allow the ability to use several types
 * of reactions without the ODEManager having to worry about how to calculate
 * rates of change differently in each case. It owns the calculation of rates,
 * which the ODEManager merely calls on it to perform.
 *
 * Each reaction holds a ReactionType variable describing its type and telling 
 * it which operatons to perform in taking input concentrations and outputing
 * their rates of change. This enables both simple interactions (like the ones
 * used by Hakim) and complicated interactions (like Hill's functions) to be
 * seamlessly added and used in different simulations.
 *
 * Currently, the class allows interactions involving one (e.g. degradation), 
 * two (e.g. phosphorylation), or three (e.g. combination) participants.
 *
 */

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

# include "old/numeric/dvec.h"
# include "ReactionType.h"
# include "GenomeReaction.h"
# include "Genome.h"

#ifndef ODEREACTION_H
#define ODEREACTION_H

class ODEReaction {

public:
	ODEReaction();
	ODEReaction( Genome& genome , int iReac , int iCell );
	ODEReaction( ReactionType type , 
				int numPart , 
				int cellPartZero , int cellPartOne , int cellPartTwo ,
				int iPartZero , int iPartOne , int iPartTwo ,
				int dxdtPartZero , int dxdtPartOne , int dxdtPartTwo , 
				double kineticZero , double kineticOne );
	~ODEReaction();
	
	int getNumPart();
	int getICell(int partNum);
	int getIPart(int partNum);
	double getDxDt(int partNum);
	
	void react( const std::vector< dvec* >& currTissue);
	
private:
	
	/* For reactions with less than three participants, all indexes and 
	 * locations for out-of-scope participants will be set to NEXIST.
	 */
	static const int NEXIST = -1;
	
	/* ODEReaction will calculate rates differently depending on its type */
	ReactionType type;
	
	/* Number of Participants: must be betwee 1 and 3 */
	int numPart;
	
	/* Which cell in the tissue is this participant in. */
	/* By allowing different cells for different participants, we enable
	 * trans-membrane reactions.
	 * The value can also be NEXIST.
	 */
	int cellPartZero;
	int cellPartOne;
	int cellPartTwo;
	
	/* Where in the ODEManager are the participants located */
	/* The ODEManager holds concentrations in a dvec. This tells us the
	 * location of each participant. With the cellPart____ variables
	 * we know the exact location of each participant
	 */
	int iPartZero;
	int iPartOne;
	int iPartTwo;
	
	/* What is the current rate of change of concentrations for the participants
	 * in this reaction */
	double dxdtPartZero;
	double dxdtPartOne;
	double dxdtPartTwo;
	
	/* Kinetic constants for the reaction. Depending on the type of reaction,
	 * these may have different interpretations.
	 */
	double kineticZero;
	double kineticOne;
	
	/* Copy and assignment operators made private */
	ODEReaction( ODEReaction *newOne ) {}
	ODEReaction& operator=( const ODEReaction& rhs ) { return *this; }
	
};

#endif