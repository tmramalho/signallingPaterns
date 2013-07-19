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
 * Default implementation of constructor for use during development.
 *
 * We construct a one cell tissue with three substances of concentration 3.0
 * each, with the default ODEReaction describing their interaction.
 *
 * The default ODEReaction is:
 *		type: COMBINATION
 *		moleculeZero + moleculeOne <--> moleculeTwo
 *		forwardRate = 0.7 * [moleculeZero] * [moleculeOne]
 *		backwardRate = 0.3 * [moleculeTwo]
 *
 */
ODEManager::ODEManager() {
	
	/* Set initial conditions and fill iTissue and currTissue with them */
	dvec* insertZero = new dvec(3,3.0);
	dvec* insertOne = new dvec(3,3.0);
	dvec* insertTwo = new dvec(3,3.0);
	
	iTissue.push_back(insertZero);
	currTissue.push_back(insertOne);
	dxdt.push_back(insertTwo);
	
	/* Initialize time */
	time = 0;
	
	/* Add the default reaction to our reaction vector. */
	reactions.push_back(new ODEReaction());
	
	/* Set dxdt vector according to reactions */
	updateRates();
}

/* Constructor: ODEManager(genome)
 * -------------------------------------------------------------------------- 
 * Constructor that reads the the information contained in a genome and uses
 * it to construct the corresponding ODE system representing it.
 */
ODEManager::ODEManager(Genome& genome) {
	
	int numMol = genome.getNumMol();
	iTissue.push_back( new dvec(numMol,3.0) );
	currTissue.push_back( new dvec(numMol,3.0) );
	dxdt.push_back( new dvec(numMol,3.0) );
	
	time = 0.0;
	
	readInReactions(genome);
	
	updateRates();
	
}


/* Destructor: ODEManager() 
 * -------------------------------------------------------------------------- 
 */
ODEManager::~ODEManager() {
	
	for ( int i = 0 ; i < iTissue.size() ; i++ ) {
		delete iTissue.at(i);
		delete currTissue.at(i);
		delete dxdt.at(i);
	}
	
	for ( int i = 0 ; i < reactions.size() ; i++ ) {
		delete reactions.at(i);
	}
	
}

/* Public Method: run(mode)
 * -------------------------------------------------------------------------- 
 * Runs the ODE by filling in the currTissue vector based on iTissue, 
 * reactions, and mode.
 */
void ODEManager::run(string mode) {
	if (mode == "rk1_det_ti"){
		rk1_det_ti ( 100 , .1 );
	}
}

/* Private Method: updateRates()
 * -------------------------------------------------------------------------- 
 * Updates the dxdt vector to contain the deterministic rate of change of all 
 * proteins given the current reactions and concentrations of proteins.
 */
void ODEManager::updateRates() {
	
	/* Reset all rates to zero */
	for ( int iCell = 0 ; iCell < currTissue.size() ; iCell++ ) {
		dxdt.at(iCell)->zero();
	}
	
	/* Update dxdt vector with new rates, which the ODEReaction objects calculates
	 * for us.
	 */
	for ( unsigned int iReac = 0 ; iReac < reactions.size() ; iReac++ ) {
		
		reactions.at(iReac)->react(currTissue);
		
		/* For this reaction, we go through the participants and update our dxdt
		 * dvecs accordingly.
		 */
		ODEReaction *currReac = reactions.at(iReac);
		
		currReac->react(currTissue); /* This is where the actual calculation occurs.
									  * Everything else is simply looking up.
									  */
		
		for ( int p = 0 ; p < reactions.at(iReac)->getNumPart() ; p++ ) {
			dxdt.at(currReac->getICell(p))->at(currReac->getIPart(p)) += currReac->getDxDt(p);
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
	for ( unsigned int iCell = 0 ; iCell < currTissue.size() ; iCell++ ) {
		dvec* cellVec = currTissue.at(iCell);
		for ( unsigned int i = 0; i < cellVec->size(); i++ ) {
			double update = (dxdt.at(iCell)->at(i))*dt;
			cellVec->at(i) += update;
		}
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
		std::cout << step << ": " << currTissue.at(0)->at(0)
		<< " " << currTissue.at(0)->at(1) << " " <<
		currTissue.at(0)->at(2) << std::endl;
		rk1_det_ti_step( dt );
	}
	std::cout << std::endl;
	std::cout << "Final Configuration" << std::endl;
	std::cout << "iTissue: " << iTissue.at(0)->at(0)
	<< " " << iTissue.at(0)->at(1) << " " <<
	iTissue.at(0)->at(2) << std::endl;
	std::cout << "currTissue: " << currTissue.at(0)->at(0)
	<< " " << currTissue.at(0)->at(1) << " " <<
	currTissue.at(0)->at(2) << std::endl;
	std::cout << "dxdt: " << dxdt.at(0)->at(0)
	<< " " << dxdt.at(0)->at(1) << " " <<
	dxdt.at(0)->at(2) << std::endl;
}

/* Private Method: readInReactions(genome)
 * -------------------------------------------------------------------------- 
 * Helper method for translation of genome.
 * 
 * We iterate through the GenomeReactions in our genome, and for each one we
 * create a corresponding ODEReaction to add to our reactions vector for use
 * during ODEIntegration.
 *
 * A main difference between ODEReactions and GenomeReactions is that while
 * our genome distinguishes between genes and proteins, our ODE does not. 
 * To our ODE, everything is simply given a slot in our dvec. 
 * 
 * Thus, in order to translate from the genome to the ODE system, we need a 
 * protocol for how we arrange our substances in our dvecs. This protocol is
 * contained in the operation of the Genome class, which can take it's
 * reactantions and tell us where each participant should be in the dvecs of
 * our ODEManager. Because the genome owns this placement operation, it will
 * be consistant everywhere in our code.
 */
void ODEManager::readInReactions(Genome& genome) {
	
	for ( int iCell = 0 ; iCell < currTissue.size() ; iCell++ ) {
		for ( int iReac = 0 ; iReac < genome.getNumReac() ; iReac++ ) {
			reactions.push_back(new ODEReaction(genome,iReac,iCell));
		}
	}	
	
}












