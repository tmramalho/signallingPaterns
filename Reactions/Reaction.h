/*
 *  Reaction.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/18/13.
 *
 */

/* Reaction class
 * -------------------------------------------------------------------------- 
 * The Reaction class is an abstract class.
 *
 * All reactions have a number of participants, which is the number of
 * molecules whose concentration gets changed due to the reaction.
 *
 * CAUTION: This is different than the number of molecules whose 
 * concentrations play a role in the rate of the reaction.
 *
 */

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
# include <vector>

# include "../old/numerics/dvec.h"
# include "ReactionType.h"
# include "../IntegrationType.h"

#ifndef REACTION_H
#define REACTION_H

class Reaction {

public:
	Reaction();
	~Reaction();
	
	int getNumPart();
	virtual int getIPart(int partNum) = 0;
	virtual double getDx(int partNum) = 0;
	
	virtual void react( std::vector< dvec* >& currTissue , std::vector< std::vector<int>* >& neighbors,
					   int iCurrCell , IntegrationType mode , double dt ) = 0;
	
protected:
	
	/* For reactions with less than three participants, all indexes and 
	 * locations for out-of-scope participants will be set to NEXIST.
	 */
	static const int NEXIST = -1;
	
	/* Reaction will calculate rates differently depending on its type */
	ReactionType type;
	
	int numPart;
	
private:
	
	/* Copy and assignment operators made private */
	Reaction( Reaction *newOne ) {}
	Reaction& operator=( const Reaction& rhs ) { return *this; }
	
};

#endif