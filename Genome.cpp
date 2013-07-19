/*
 *  Genome.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/17/13.
 *
 */

#include "Genome.h"


/* Constructor: Genome()
 * -------------------------------------------------------------------------- 
 * Default implementation of constructor for use during development.
 *
 * A genome consisting of three proteins and one reaction between them. This
 * reaction is:
 *
 *		type: COMBINATION
 *		proteinZero + proteinOne <--> proteinTwo
 *		forwardRate = 0.7 * [proteinZero] * [proteinOne]
 *		backwardRate = 0.3 * [proteinTwo]
 *
 */

Genome::Genome() {
	reactions.push_back(new GenomeReaction());
	
	proteins.push_back(new Protein(0,NEXIST,NEXIST));
	proteins.push_back(new Protein(1,NEXIST,NEXIST));
	proteins.push_back(new Protein(2,NEXIST,NEXIST));
	
	proteins.at(0)->addReaction(reactions.at(0));
	proteins.at(1)->addReaction(reactions.at(0));
	proteins.at(2)->addReaction(reactions.at(0));
}

/* Destructor: ~Genome()
 * -------------------------------------------------------------------------- 
 */
Genome::~Genome() {
	for ( int i = 0 ; i < reactions.size() ; i++ ) {
		delete reactions.at(i);
	}
	for ( int i = 0 ; i < genes.size() ; i++ ) {
		delete genes.at(i);
	}
	for ( int i = 0 ; i < proteins.size() ; i++ ) {
		delete proteins.at(i);
	}
}

/* Public Method: getNumMol()
 * -------------------------------------------------------------------------- 
 * Returns the number of molecules, including genes, proteins, and complexes,
 * in this genome.
 */
int Genome::getNumMol() {
	return genes.size() + proteins.size();
}

/* Public Method: getNumReac()
 * -------------------------------------------------------------------------- 
 * Returns the number of reactions in this genome.
 */
int Genome::getNumReac() {
	return reactions.size();
}

/* Public Method: getReacRef(i)
 * -------------------------------------------------------------------------- 
 * Returns a reference to the reaction held in this genome with index i.
 */
GenomeReaction* Genome::getReacRef(int i) {
	return reactions.at(i);
}

/* Public Method: getIPartForODE( int reacNum , int partNum )
 * -------------------------------------------------------------------------- 
 * Returns the index where the partNum participant of the reacNum reaction
 * should be stored in the ODEManager.
 *
 * The current convention is to first list all of the genes in the genome in
 * the order in which they appear in the genes vector, followed by all the 
 * proteins in the order in which they appear in that vector.
 *
 * See "Note on translation to ODEManager" in the description of the Genome
 * class in Genome.h to understand why this method is necessary.
 *
 */
int Genome::getIPartForODE( int iReac , int partNum ) {
// Make sure that for particles that don't exist, you return index NEXIST.	
	GenomeReaction* reacRef = reactions.at(iReac);
	ReactionType reacType = reacRef->getReacType();
	
	switch (reacType) {
		case PROMOTION:
			break;
		case DEGRADATION:
			break;
		case COMBINATION:
			return genes.size() + reactions.at(iReac)->getIPart(partNum);
			break;
	}	
}






