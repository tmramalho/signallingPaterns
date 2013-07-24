/*
 *  Molecule.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/17/13.
 *
 */

/* Molecule class
 * -------------------------------------------------------------------------- 
 * The Molecule class is a general class from which Gene and Protein inherit.
 * Molecules are held in the genome as members of a gene network. We want 
 * the molecules to be aware of where they are in the genome, so we give each
 * of them knowledge of their index in the vector they are held in in the 
 * genome and a reference to the reactions they take part in.
 *
 */


#ifndef MOLECULE_H
#define MOLECULE_H

# include <vector>
# include "../Reactions/Reaction.h"

using namespace std;

class Molecule {
	
public:
	Molecule();
	~Molecule();
	
	unsigned int getISelf();
	double getInitConc();
	void addReaction(Reaction* const reacRef);
	
	virtual void updateIndices( int firstIndex , int numInsertion) = 0;
	
protected:
	
	static const int NEXIST = -1;
	
	unsigned int iSelf;
	
	double initConc;
	
	/* The GenomeReactions are owned by the Genome class. The Molecule
	 * class references them only so the molecule can know which
	 * reactions it participates in without having to search for itself
	 * in all reactions. 
	 *
	 * The molecule destructor should NOT delete the GenomeReactions.
	 *
	 */
	vector<Reaction*> reactions;
	
};

#endif