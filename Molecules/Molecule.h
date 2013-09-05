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
# include <cstring>
# include <cstdlib>
# include <iostream>
# include <fstream>
# include <algorithm>
# include <set>
# include <iterator>
# include "../Reactions/Reaction.h"
# include "../old/helpers/SettingsCont.h"
# include "../nexist.cpp"

using namespace std;

class Molecule {
	
public:
	Molecule( const Molecule& other );
	~Molecule() {}
	
	int get_i_self();
	double get_init_conc();
	void add_reaction( int i_reac );
	
	virtual void update_mol_indices( int first_index , int num_insertion ) = 0;
	template <class UpdateClass>
	void update_reac_indices( UpdateClass update_index );
	
	virtual void print_info ( std::string line_start ) = 0;
	virtual void to_file (ofstream& file,std::string line_start) = 0;
	
	/* Iterator through reactions molecule is relevant to */
	std::set<int>::iterator reacs_begin() const;
	std::set<int>::iterator reacs_end() const;
	
	
	
protected:
	
	Molecule();
	
	SettingsCont* _sc_ref;
	
	int _i_self;
	
	double _init_conc;
	
	/* The Reactions referenced in the _reactions vector are owned by the 
	 * Manager class. The Molecule class references them only so the molecule
	 * can know which reactions it participates in without having to search 
	 * for itself in all reactions. 
	 *
	 * The molecule destructor should NOT delete the Reactions.
	 *
	 */
	std::set<int> _reactions;
	
private:
	Molecule& operator=( const Molecule& rhs );
	Molecule(Molecule* newOne) {}
};

#endif