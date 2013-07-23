/*
 *  Manager.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/12/13.
 *
 */

#include "Manager.h"

/* Constructor: Manager()
 * -------------------------------------------------------------------------- 
 * Default implementation of constructor for use during development.
 *
 * Currently attempts to create a Delta-Notch like network.
 */

Manager::Manager() {

	/* Create neighbors vector: Linear tissue three long */
	/* < [2] , [1,3] , [2] > */
	neighbors.push_back(new vector<int>(1,1));
	vector<int>* addition = new vector<int>;
	addition->push_back(0);
	addition->push_back(2);
	neighbors.push_back(addition);
	addition = new vector<int>(1,1);
	neighbors.push_back(addition);
	
	/* Initialize tissue vectors */
	/* The substances in our genome are, in order, 
	 * 
	 *		[ d   n   d:N   d:N:N   D   N]
	 *
	 * Where lower-case refers to a gene, upper-case to a protein.
	 *
	 * Initial concentrations are
	 *
	 *		Cell 1: [ 1.0  1.0  0.0  0.0  0.0  1.0 ]
	 *		Cell 2: [ 1.0  1.0  0.0  0.0  1.0  0.0 ]
	 *		Cell 3: [ 1.0  1.0  0.0  0.0  0.0  1.0 ]
	 *
	 */
	double nums[5];
	nums[0] = 1;
	nums[1] = 1;
	nums[2] = 0;
	nums[3] = 0;
	nums[4] = 0.0;
	nums[5] = 1.0;
	iTissue.push_back(new dvec(6,nums));
	nums[4] = 1.0;
	nums[5] = 0.0;
	iTissue.push_back(new dvec(6,nums));
	nums[4] = 0.0;
	nums[5] = 1.0;
	iTissue.push_back(new dvec(6,nums));
	
	for ( int i = 0 ; i < 3 ; i++ ) {
		currTissue.push_back(new dvec(6,iTissue.at(i)));
		dx.push_back(new dvec(6,0.0));
	}
	
	/* Initialize gene and protein vectors */
	/* Reference detailed constructors and array of
	 * elements of genome above to understand 
	 * paramaters.
	 */
	genes.push_back(new Gene(0,4,-1,-1));
	genes.push_back(new Gene(1,5,-1,-1));
	genes.push_back(new Gene(2,4,5,0));
	genes.push_back(new Gene(3,4,5,2));
	
	proteins.push_back(new Protein(4,-1,-1));
	proteins.push_back(new Protein(5,-1,-1));
	
	/* Initialize time variables */
	dt = .01;
	time = 0.0;
	
	/* Initialize Integration Mode */
	mode = RK1_DET_TI;
	
	/* Initialize Reactions */
	/* d --> D with kinetic constant 1.0 */
	reactions.push_back(new PromReaction(0,4,0.0,1.0));
	/* d:N --> D with kinetic constant 1.0 */
	reactions.push_back(new PromReaction(2,4,0.0,1.0));
	/* d:N:N --> D with kinetic constant 0 */
	reactions.push_back(new PromReaction(3,4,0.0,0));
	/* n --> N with kinetic constant 0 */
	reactions.push_back(new PromReaction(1,5,0.0,0.0));
	/* d + N <--> d:N with forwardKinetic 1, backwardKinetic 1 */
	reactions.push_back(new PromBindingReaction(0,2,4,
												0.0,0.0,
												1.0,1.0));
	/* d:N + N <--> with forwardKinetic 1, backwardKinetic 1 */
	reactions.push_back(new PromBindingReaction(2,3,4,
												0.0,0.0,
												1.0,1.0));
	/* D promotes N, with kinetic = 2.0, K = .35 */
	reactions.push_back(new LatPromReaction(5,4,
											0.0,2.0,
											.35));
	/* D --> nothing with kinetic constant 1.0 */
	reactions.push_back(new DegReaction(4,0.0,1.0));
	/* N --> nothing with kinetic constant 2.0 */
	reactions.push_back(new DegReaction(5,0.0,2.0));
	
	updateDx();
	
}

/* Destructor: ~Manager() 
 * -------------------------------------------------------------------------- 
 */
Manager::~Manager() {
	
	for ( int i = 0 ; i < iTissue.size() ; i++ ) {
		delete iTissue.at(i);
		delete currTissue.at(i);
		delete dx.at(i);
	}
	
	for ( int i = 0 ; i < reactions.size() ; i++ ) {
		delete reactions.at(i);
	}
	
	for ( int i = 0 ; i < genes.size() ; i++ ) {
		delete genes.at(i);
	}
	
	for ( int i = 0 ; i < proteins.size() ; i++ ) {
		delete proteins.at(i);
	}
	
	neighbors.erase(neighbors.begin(),neighbors.begin()+2);
	
}

/* Public Method: integrate(mode,dt,numSteps)
 * -------------------------------------------------------------------------- 
 * Runs the ODE by filling in the currTissue vector based on iTissue, 
 * reactions, and mode.
 */
void Manager::integrate( int numSteps ) {
	for ( int step = 0 ; step < numSteps ; step++ ) {
		
		/* For the User */
		printState();
		
		/* Update currTissue according to current rates */
		/* Our dx vector always should have the appropriate update for the 
		 * next step.
		 */
		for ( unsigned int iCell = 0 ; iCell < currTissue.size() ; iCell++ ) {
			dvec* currCell = currTissue.at(iCell);
			for ( unsigned int iMol = 0; iMol < currCell->size(); iMol++ ) {
				double update = dx.at(iCell)->at(iMol);
				currCell->at(iMol) += update;
			}
		}
		
		/* Never return a public method without dx vector being up to date */
		updateDx();
		
		/* Progress time */
		time += dt;
	}
}

/* Private Method: updateDx()
 * -------------------------------------------------------------------------- 
 */
void Manager::updateDx() {
	
	/* Reset all rates to zero */
	for ( int iCell = 0 ; iCell < currTissue.size() ; iCell++ ) {
		dx.at(iCell)->zero();
	}
	
	for (unsigned int iCell = 0; iCell < iTissue.size() ; iCell++) {
		/* Update dx vector with new rates, which the Reaction objects calculates
		 * for us.
		 */
		for ( unsigned int iReac = 0 ; iReac < reactions.size() ; iReac++ ) {
		
			/* For this reaction, we go through the participants and update our dxdt
			 * dvecs accordingly.
			 */
			Reaction *currReac = reactions.at(iReac);
			
			/* Have reaction object calculate change */
			currReac->react(currTissue,neighbors,iCell,mode,dt); 
			
			/* Transfer changes to dx vector */
			for ( int partNum = 0 ; partNum < reactions.at(iReac)->getNumPart() ; partNum++ ) {
				double update = currReac->getDx(partNum);
				dx.at(iCell)->at(currReac->getIPart(partNum)) += update;
			}
		}
	}
}

/* Public Method: updateDx()
 * -------------------------------------------------------------------------- 
 * For the user.
 */
void Manager::printState() {
	std::cout << "time: " << time << std::endl;
	for ( int i = 0 ; i < currTissue.size() ; i++ ) {
		std::cout << "Cell: " << i << std::endl;
		std::cout << currTissue.at(i)->at(0) << " :::: " <<
		currTissue.at(i)->at(1) << " :::: " << currTissue.at(i)->at(2) <<
		" :::: " << currTissue.at(i)->at(3) << " :::: " << currTissue.at(i)->at(4) <<
		" :::: " << currTissue.at(i)->at(5) << std::endl;
	}
	std::cout << std::endl;
}





