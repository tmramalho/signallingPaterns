/*
 *  ODEManager.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/12/13.
 *
 */

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

# include "old/numeric/dvec.h"
# include "ODEReaction.h"

# ifndef _ODEManager_H
# define _ODEManager_H

using namespace std;

class ODEManager {
	
public:
	
	ODEManager();
	~ODEManager();
	
	void run(string mode);
	
private:
	
	std::vector< dvec* > iTissue, currTissue, dxdt;
	/* At all times, it should have the current rates of
	 * change of concentrations of proteins. Implementation
	 * of all methods should be such to maintain this.
	 */ 
	double time;
	
	/* Vector of ODEReactions */
	std::vector< ODEReaction* > reactions;
	
	/* Copy and assignment operators made private */
	ODEManager( ODEManager *newOne );
	ODEManager& operator=( const ODEManager& rhs );
	
	/* Helpers for running of ODEs */
	void updateRates();
	void rk1_det_ti_step ( double dt );
	void rk1_det_ti ( int numSteps , double dt );
	
};


# endif