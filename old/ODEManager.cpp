/*
 *  ODEManager.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/12/13.
 *
 */

#include "ODEManager.h"

/* Constructor: ODEManager()
 * -------------------------------------------------------------------------- 
 * Default implementation of constructor for use during development 
 */
ODEManager::ODEManager():iTissue(1),currTissue(1),dxdt(1) {
	
	/* Set initial conditions and fill iTissue and currTissue with them */
	dvec insertOne(3,3.0);
	iTissue.at(0).set(insertOne);
	currTissue.at(0).set(insertOne);
	dxdt.at(0).set(insertOne); // We need to make sure the dxdt begins with a dvec of the appropriate size.
	
	/* Initialize time */
	time = 0;
	
	/* Create reaction linked list */
	/* Currently we initialize to a list of length one with the parameters as 
	 * listed.
	 */
	reactions = new Reaction[1];
	reactions[1].reactantOne = 0;
	reactions[1].reactantTwo = 1;
	reactions[1].product = 2;
	reactions[1].forwardRate = .3;
	reactions[1].backwardRate = .7;
	
	/* Set dxdt vector according to reactions */
	updateRates();
}

/* Constructor: ODEManager() 
 * -------------------------------------------------------------------------- 
 */
ODEManager::~ODEManager() {
	
	/* Delete reactions array */
	delete[] reactions;
	
}

/* Copy Operator
 * -------------------------------------------------------------------------- 
 */
ODEManager::ODEManager( ODEManager *newOne ) {
	
}

/* Assignment Operator 
 * -------------------------------------------------------------------------- 
 */
ODEManager& ODEManager::operator=( const ODEManager& rhs ) { return *this; }

/* Public Method: run(mode)
 * -------------------------------------------------------------------------- 
 * Runs the ODE by filling in the fTissue vector based on iTissue, Reactions,
 * and mode.
 */
/* Temp: mode is irrelevant, but will specify which numerical method we are
 * going to use. For the moment, we use rk1_det_ti.
 */
void ODEManager::run(string mode) {
	rk1_det_ti ( 100 , .1 );
}

/* Private Method: updateRates()
 * -------------------------------------------------------------------------- 
 * Updates the dxdt vector to contain the deterministic rate of change of all 
 * proteins given the current reactions and concentrations of proteins.
 */
void ODEManager::updateRates() {
	
	for ( int iCell = 0 ; iCell < currTissue.size() ; iCell++ ) {
		/* Reset all rates to zero */
		dxdt.at(iCell).zero();
		
		/* Find dvec of protein concentrations in this cell. Then iterate 
		 * through the reactions to update rates of change of these 
		 * concentrations
		 */
		for ( int iReaction = 0; iReaction < 1 ; iReaction++ ) { // for now we have only one cell
			
			/* forwardFlow is lambdaOne * [reactantOne] * [reactantTwo] */
			double forwardFlow = reactions[iReaction].forwardRate * 
			currTissue.at(iCell).at(reactions[iReaction].reactantOne) * 
			currTissue.at(iCell).at(reactions[iReaction].reactantTwo);
			/* backwardFlow is lambdaTwo * [product] */
			double backwardFlow = reactions[iReaction].backwardRate *
			currTissue.at(iCell).at(reactions[iReaction].product);
			
			/* Update concentrations in our dvec */
			dxdt.at(iCell).at(reactions[iReaction].reactantOne) += backwardFlow - forwardFlow;
			dxdt.at(iCell).at(reactions[iReaction].reactantTwo) += backwardFlow - forwardFlow;
			dxdt.at(iCell).at(reactions[iReaction].product) += forwardFlow- backwardFlow;
			
		}
	}	
}

/* Private Method: rk1_det_ti_step(dt)
 * -------------------------------------------------------------------------- 
 * Updates currState by one time step dt, according to deterministic, first-
 * order Runge-Kutta.
 */
void ODEManager::rk1_det_ti_step ( double dt ) {
	
	/* Update currTissue according to current rates */
	for ( int iCell = 0 ; iCell < currTissue.size() ; iCell++ ) {
		currTissue.at(iCell) += (dxdt.at(iCell))*dt;
	}
	
	/* Update dxdt vector to maintain accuracy of current representation*/
	updateRates();
	
	/* Update time */
	time += dt;
}

/* Private Method: rk1_det_ti(numSteps,dt)
 * -------------------------------------------------------------------------- 
 * Updates currState by numSteps time steps dt, according to deterministic 
 * first-order Runge-Kutta.
 */
void ODEManager::rk1_det_ti ( int numSteps , double dt ) {
	for ( int step = 0 ; step < numSteps ; step++ ) {
		std::cout << step << ": " << currTissue.at(0).at(0) << std::endl;
		rk1_det_ti_step( dt );
	}
	std::cout << currTissue.at(1).at(0) << " " <<
	currTissue.at(1).at(1) << " " <<
	currTissue.at(1).at(2) << " " << std::endl;
	
}



