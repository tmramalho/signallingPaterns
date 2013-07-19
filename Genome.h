/*
 *  Genome.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/17/13.
 *
 */

/* Genome class
 * -------------------------------------------------------------------------- 
 * The Genome class contains the detailed information about our gene network
 * in a format that is easy to mutate. (This functionality will be added
 * later.)
 *
 * The genome class is basically a container of all the Genes, Proteins, and
 * GenomeReactions that constitute a gene network. These are themselves
 * classes which contain the relevant information to their structure. See
 * the classes individually for what this consists of.
 *
 * The Genome can access these participants based on type and provide info
 * about them.
 * 
 * Note on translation to ODEManager:
 * Our genome consists of a certain number of genes and proteins in a gene
 * network. When an ODEManager is constructed, these molecules must be stored
 * in the ODEManager in some consistent order. The genome controls what this
 * order is, and for any substance can be asked its index in the list of 
 * molecules. The ODEManager should not determine this order independently
 * to avoid inadvertant inconsistancies. The order is only made available to
 * the ODEManager by asking the genome for it.
 *
 */

#ifndef GENOME_H
#define GENOME_H

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
# include <vector>

# include "Gene.h"
# include "Protein.h"
# include "Molecule.h"
# include "GenomeReaction.h"

class Genome {
	
public:
	
	Genome();
	~Genome();
	
	int getNumMol();
	int getNumReac();	
	
	GenomeReaction* getReacRef(int i);
	
	/* Ordering functionality */
	/* Tells the index in the ODEManager the the partNum participant in the reacNum
	 * reaction should be in.
	 */
	int getIPartForODE( int iReac , int partNum );
	
private:

	static const int NEXIST = -1;
	
	/* Vector of reactions pointers */
	std::vector< GenomeReaction* > reactions;

	std::vector<Gene*> genes;
	std::vector<Protein*> proteins;

};

#endif



