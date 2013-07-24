/*
 *  Gene.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/17/13.
 *
 */

/* Gene class, subclass of Molecule
 * -------------------------------------------------------------------------- 
 * In addition to the information a Gene objects inherits from the Molecule
 * class, the gene has knowledge of which protein it codes for. 
 *
 * A gene object can also represent a gene-promoter complex, in which case it
 * has knowledge of which protein is bound to it's promoter and which gene it 
 * is when not in a complex (root). This information is stored only via the
 * index of these proteins or genes in the genome data structure.
 *
 */

# include "Molecule.h"

#ifndef GENE_H
#define GENE_H

class Gene : public Molecule {

public:
	Gene( int iSelf, int iProduct , int iBoundPromoter , int iRoot , double initConc );
	~Gene();
	
	int getIProduct();
	int getIBoundPromoter();
	int getIRoot();
	
	virtual void updateIndices( int firstIndex , int numInsertions );
	
protected:
	
	/* All indices refer to their index in the genome containing them. */
	
	/* Which protein does the gene code for */
	int iProduct;
	
	/* If the gene has a promoter bound to it, which is it */
	int iBoundPromoter;
	
	/* If the gene has a promoter bound to it, which is the base, promoter free
	 * gene.
	 */
	int iRoot;
	
};

#endif