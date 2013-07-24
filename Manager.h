/*
 *  Manager.h
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
# include "Reactions/Reaction.h"
# include "Reactions/CombReaction.h"
# include "Reactions/DegReaction.h"
# include "Reactions/PromBindingReaction.h"
# include "Reactions/PromReaction.h"
# include "Reactions/LatPromReaction.h"
# include "Reactions/PromBindingReaction.h"
# include "IntegrationType.h"
# include "Molecules/Molecule.h"
# include "Molecules/Protein.h"
# include "Molecules/Gene.h"

# ifndef _Manager_H
# define _Manager_H

using namespace std;

/* Manager class
 * -------------------------------------------------------------------------- 
 * This class holds all the necessary information for a particular genome 
 * and the integration of the corresponding ODE.
 *
 * It also constantly stores information about which integration mode it
 * it should use to integrate, the current time in the integration, and the
 * time step it should use to integrate.
 *
 * This means that we just call integrate() for the integration, and must
 * change these settings before calling integrate() if we want to manage
 * them.
 *
 */

class Manager {
	
public:
	
	Manager();
	~Manager();
	
	void setMode( IntegrationType mode );
	IntegrationType getMode();
	
	void setDt( double dt );
	double getDt();
	
	/* Adding Reactions */
	void addGene( double prodRate , double degRate );
	void addPromBindingReac( int iGeneInGenes , int iProtInProts ,
							double forwardKinetic , double backwardKinetic , 
							double newProdRate );
	//void addPhosphReac();
	//void addPartialDegReac();
	void addCombReaction( int iProtZeroInProts , int iProtOneInProts , 
						 double forwardKinetic , double backwardKinetic );
	//void addPartialCatDegReac();
	void addLatPromReac( int iLocProtInProts , int iNeighborProtInProts ,
						double kinetic , double K );
	
	/* Functions for running the ODE */
	void initialize();
	void integrate( int numSteps );
	
	/* Display Function */
	void printState();
	//void printGenome();
	
private:
	
	static const int NEXIST = -1;
	
	/* Neighbors info */
	/* Our tissue data is held in a vector, with each element a separate
	 * cell. Because our genetic network will include membrane reactions,
	 * we would like each cell to be able to know who its neighbors are.
	 * Each element of this vector holds a vector of the indices of the
	 * corresponding cells.
	 *
	 * For example, if we have a 2x2 tissue, we might index the cells
	 * from left to right, top to bottom, as
	 *
	 *		1	2
	 *
	 *		3	4
	 *
	 * Then, our vector of neighbors would look like this:
	 *
	 *		< [2,3] , [1,4] , [1,4] , [2,3] >
	 *
	 * This is how we store the geometry of our lattice, and it can 
	 * accomodate any geometry we desire.
	 * Any use can then easily use the manager class to study any geometry
	 * he or she wants.
	 *
	 */
	std::vector< std::vector<int>* > neighbors;
	
	/* State vectors */
	/* The class is implemented so that these always up-to-date. In particular,
	 * no public method should return before the dxdt vector has the correct
	 * rates given the state of currTissue. Then, to get the current rates of
	 * change in our ODE, we only need to ask the ODEManager to look in the
	 * dxdt vector. All calculations should have already been done.
	 */
	std::vector< dvec* > iTissue, currTissue, dx;
	
	/* Gene and protein info */
	std::vector<Gene*> genes;
	std::vector<Protein*> proteins;
	
	/* Time variables */
	double dt; 
	double time;
	
	/* Number of molecules */
	int numMol;
	
	/* Vector of Reactions */
	std::vector< Reaction* > reactions;
	
	/* Integration Mode */
	IntegrationType mode;
	
	/* Copy and assignment operators made private */
	Manager( Manager *newOne );
	Manager& operator=( const Manager& rhs );
	
	/* Helpers for running of ODEs */
	void updateDx();
	
	/* Helpers for mutating genome */
	void resizeDVecs();
	void updateIndices( int firstIndex , int numInsertions );
	
};


# endif