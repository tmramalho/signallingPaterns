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

#ifndef PROTEIN_H
#define PROTEIN_H

# include "Molecule.h"

class Protein : public Molecule {

public:
	Protein( int i_self , int i_root_zero , int i_root_one , double init_conc );
	Protein( const Protein& other );
	Protein( std::ifstream& file );
	~Protein() {}
	
	int get_i_root_zero();
	int get_i_root_one();
	
	void update_mol_indices( int first_index , int num_insertion );

	void print_info ( std::string line_start );
	void to_file ( std::ofstream& file , std::string line_start);
	
protected:
	
	/* All indices refer to their index in the genome containing them. */
	
	/* If the protein is a complex, what are the indices of the proteins that
	 * make it up. (Note these can be complexes themselves, or be the same.)
	 */
	int _i_root_zero;
	int _i_root_one;
	
	/* File Format Variables */
	static const int NUM_LINE = 4;
	static const int I_SELF_LINE = 0;
	static const int I_ROOT_ZERO_LINE = 1;
	static const int I_ROOT_ONE_LINE = 2;
	static const int INIT_CONC_LINE = 3;
	
private:
	Protein(Protein* newOne) {}
	Protein& operator=( const Protein& rhs ) {return *this;}

};

#endif