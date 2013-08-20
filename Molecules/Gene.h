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
	Gene( int i_self, int i_product , int i_bound_promoter , int i_root , double init_conc );
	~Gene() { std::cout << "Gene Default Destructor\n"; }
	
	int get_i_product();
	int get_i_bound_promoter();
	int get_i_root();

	virtual void update_indices( int first_index , int num_insertions );
	virtual void print_info ( std::string line_start );
	virtual void to_file ( std::ofstream& file , std::string line_start);
	
protected:
	
	/* All indices refer to their index in the genome containing them. */
	
	/* Which protein does the gene code for */
	int _i_product;
	
	/* If the gene has a promoter bound to it, which is it */
	int _i_bound_promoter;
	
	/* If the gene has a promoter bound to it, which is the base, promoter free
	 * gene.
	 */
	int _i_root;
	
private:
	Gene(Gene* newOne) {}
	Gene& operator=( const Gene& rhs ) {return *this;}
	
};

#endif