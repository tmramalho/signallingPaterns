/*
 *  ODEManager.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/12/13.
 *
 */

/* ODEManager class
 * -------------------------------------------------------------------------- 
 * The ODEManager class efficiently integrates the ODE's describing our gene
 * network. It holds the following data:
 *
 *		- The initial concentration of all substances in the cells (genes,
 *		proteins, protein complexes, gene-promoter complexes).
 *		- The reactions describing the rates of change of these concentrations.
 *		- The current time in our simulation
 *		- The current concentration of all substances in the cells.
 *		- The current rate of change of all substances in the cells.
 *
 * The ODEManager has the capability to integrate using the following modes:
 *		- first order, deterministic, Runge-Kutta ("rk1_det_ti")
 *
 * It does not contain the following information not relevant to the efficient
 * integration of the ODE's, which is instead contained by the Genome class:
 *
 *		- The identify of each substance (e.g. whether it is a protein or a 
 *		gene).
 *		- The relationship between different substances. E.g. if a protein is
 *		a complex of other proteins, which proteins it is a complex of.
 *
 * The Genome class contains this extra information. The ODEManager reads from
 * a Genome object in its construction to build the corresponding ODE system.
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
# include "Genome.h"
# include "GenomeReaction.h"

# ifndef _ODEManager_H
# define _ODEManager_H

using namespace std;

class ODEManager {
	
public:
	
	ODEManager();
	ODEManager(Genome& genome);
	~ODEManager();
	
	void run(string mode);
	
private:
	
	/* State vectors */
	/* The class is implemented so that these always up-to-date. In particular,
	 * no public method should return before the dxdt vector has the correct
	 * rates given the state of currTissue. Then, to get the current rates of
	 * change in our ODE, we only need to ask the ODEManager to look in the
	 * dxdt vector. All calculations should have already been done.
	 */
	std::vector< dvec* > iTissue, currTissue, dxdt;
	
	/* Current time */
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
	
	/* Helpers for translation of genome */
	void readInReactions(Genome& genome);
};


# endif