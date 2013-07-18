/*
 *  ODEReaction.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/18/13.
 *
 */

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

# include "old/numeric/dvec.h"

enum ReactionType { PROMOTION , DEGRADATION , COMBINATION };

class ODEReaction {

public:
	ODEReaction();
	~ODEReaction();
	
	int getNumPart();
	int getCellLoc(int numParticipant);
	int getMolLoc(int numParticipant);
	double getDxDt(int numParticipant);
	
	void react(std::vector< dvec* >& currTissue);
	
private:
	
	static const int NEXIST = -1;
	
	// ODEReaction will react differently depending on the reaction type.
	ReactionType type;
	
	// There are between one and three participants in a reaction.
	int numParticipants;
	
	// We will assign NEXIST to the cell or location of any non-existant participant.
	int cellPartZero;
	int cellPartOne;
	int cellPartTwo;
	
	int locPartZero;
	int locPartOne;
	int locPartTwo;
	
	double dxdtPartZero;
	double dxdtPartOne;
	double dxdtPartTwo;
	
	double kineticZero;
	double kineticOne;
	
};