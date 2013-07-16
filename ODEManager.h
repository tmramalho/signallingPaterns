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

# include "genlib.h"
# include "dvec.h"

# ifndef _ODEManager_H
# define _ODEManager_H

class ODEManager {
	
public:
	
	ODEManager();
	~ODEManager();
	
	void run(string mode);
	
private:
	
	struct Reaction {
		int reactantOne;
		int reactantTwo;
		int product;
		double forwardRate;
		double backwardRate;
		Reaction *link;
	};
	
	std::vector< dvec* > *iTissue; /* Can this be static so we don't have to
									    * reinitalize after every call to ODE?
									    */
	std::vector< dvec* > *currTissue;
	std::vector< dvec* > *dxdt; /* At all times, it should have the current rates of
							* change of concentrations of proteins. Implementation
							* of all methods should be such to maintain this.
							*/ 
	double time;
	
	/* Linked list of reactions */
	Reaction *reactions;
	
	/* Copy and assignment operators made private */
	ODEManager( ODEManager *newOne );
	ODEManager& operator=( const ODEManager& rhs );
	
	/* Helpers for running of ODEs */
	void updateRates();
	void rk1_det_ti_step ( double dt );
	void rk1_det_ti ( int numSteps , double dt );
	
};


# endif