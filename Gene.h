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
	Gene( int iProduct , int iBoundPromoter , int iRoot );
	~Gene();
	
private:
	int iProduct;
	int iBoundPromoter;
	int iRoot;
};

#endif