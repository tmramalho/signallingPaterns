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
ODEManager::ODEManager() {
	
	/* Create iTissue vector */ // Temp: initialize to specifications
	iTissue = new std::vector< dvec* >(1);
	iTissue->at(0) = new dvec(3, 3.0);
	
	/* Initialize fTissue to all proteins at 0 concentration */
	currTissue = new std::vector< dvec* >(1);
	currTissue->at(0) = new dvec(iTissue->at(0));
	
	/* Initialize time */
	time = 0;
	
	/* Create reaction linked list */
	/* Currently we initialize to a list of length one with the parameters as 
	 * listed.
	 */
	reactions = new Reaction;
	reactions->reactantOne = 0;
	reactions->reactantTwo = 1;
	reactions->product = 2;
	reactions->forwardRate = .3;
	reactions->backwardRate = .7;
	reactions->link = NULL;
	
	/* Initialize dxdt vector. */
	dxdt = new std::vector< dvec* >(1);
	dxdt->at(0) = new dvec(3,3.0);
	updateRates();
}

/* Constructor: ODEManager() 
 * -------------------------------------------------------------------------- 
 */
ODEManager::~ODEManager() {
	
	/* Delete iTissue and what it points to. */
	for ( int i = 0; i < iTissue->size(); i++ ) {
		delete iTissue->at(i);
	}
	delete iTissue;
	
	/* Delete fTissue and what it points to. */
	for ( int i = 0; i < currTissue->size(); i++ ) {
		delete currTissue->at(i);
	}
	delete currTissue;
	
	/* Delete reactions linked list */
	Reaction *curr = reactions;
	while ( curr != NULL ) {
		Reaction *next = curr->link;
		delete curr;
		curr = next;
	}
	
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
	
	for ( int iCell = 0 ; iCell < currTissue->size() ; iCell++ ) {
		/* Reset all rates to zero */
		dxdt->at(iCell)->zero();
		
		/* Find dvec of protein concentartions in this cell. Then iterate 
		 * through the reactions to update rates of change of these 
		 * concentrations
		 */
		dvec *currCell = currTissue->at(iCell);
		Reaction *currReaction = reactions;
		while ( currReaction != NULL ) {
			
			/* forwardFlow is lambdaOne * [reactantOne] * [reactantTwo] */
			double forwardFlow = currReaction->forwardRate * 
			currCell->at(currReaction->reactantOne) * 
			currCell->at(currReaction->reactantTwo);
			/* backwardFlow is lambdaTwo * [product] */
			double backwardFlow = currReaction->backwardRate *
			currCell->at(currReaction->product);
			
			/* Update concentrations in our dvec */
			dxdt->at(iCell)->at(currReaction->reactantOne) += backwardFlow - forwardFlow;
			dxdt->at(iCell)->at(currReaction->reactantTwo) += backwardFlow - forwardFlow;
			dxdt->at(iCell)->at(currReaction->product) += forwardFlow- backwardFlow;
			
			/* Move to next reaction */
			currReaction = currReaction->link;
			
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
	for ( int iCell = 0 ; iCell < currTissue->size() ; iCell++ ) {
		(*currTissue->at(iCell)) += (*dxdt->at(iCell))*dt; // I think I might be having a problem here.
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
		std::cout << "We made it here!" << 3 << std::endl;
		rk1_det_ti_step( dt );
	}
	std::cout << currTissue->at(1)->at(0) << " " <<
	currTissue->at(1)->at(1) << " " <<
	currTissue->at(1)->at(2) << " " << std::endl;
	
}



