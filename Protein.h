/*
 *  Protein.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/17/13.
 *
 */

/* Protein class, subclass of Molecule
 * -------------------------------------------------------------------------- 
 * In addition to the information a Gene objects inherits from the Molecule
 * class, the protein, if it is a complex of proteins or the phosphorylated
 * version of a protein, contains information about which proteins it
 * derives from. This information is stored only via the index of these 
 * proteins in the genome data structure.
 *
 */

# include "Molecule.h"

#ifndef PROTEIN_H
#define PROTEIN_H

class Protein : public Molecule {

public:
	Protein( int iSelf , int iRootZero , int iRootOne );
	~Protein();
	
private:
	
	/* All indices refer to their index in the genome containing them. */
	
	/* If the protein is a complex, what are the indices of the proteins that
	 * make it up. (Note these can be complexes themselves, or be the same.)
	 */
	int iRootZero;
	int iRootOne;
};

#endif