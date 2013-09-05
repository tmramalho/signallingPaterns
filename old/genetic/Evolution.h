#ifndef _EVOLUTION_H
#define _EVOLUTION_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include "../numeric/dvec.h"
#include "../opt/FitnessFunction.h"
#include "../opt/DNScore.h"
#include "../helpers/SettingsCont.h"
#include "../../Manager.h"
#include "../../ConstructionMethod.h"
#include "omp.h"
#include "../old/genetic/EvAlgorithm.h"
#include "../old/genetic/Organism.h"
#include "../old/genetic/Lineage.h"
#include "../Molecules/Molecule.h"

# include <boost/random/uniform_int_distribution.hpp>
# include <boost/random/uniform_real_distribution.hpp>
# include <boost/random/discrete_distribution.hpp>
# include <boost/random/mersenne_twister.hpp>

/* Other Boost Libraries */
# include <boost/lexical_cast.hpp>

class Evolution {

public:
	Evolution(FitnessFunction *ff);

	virtual ~Evolution();

	void colonize( ConstructionMethod m );

	Manager *train( EvAlgorithm alg , int run_num);
	
	void genomes_to_file(int alg_num,int run_num,int gen_num);
	
	/* Data Tracking */
	static Organism* new_organism();
	void phylogeny_to_file( std::ofstream& file );
	void evaluations_to_file( std::ofstream& file );
	vector<Lineage*> create_lineages();
	
		
private:
	
	/* PRIVATE DATA MEMBERS */
	/* Population Vector */
	std::vector<Manager*> _pop;
	
	/* Algorithm Parameters */
	FitnessFunction *_ff;
	int _verbose, _ng,  _na, _mutations_per_gen;
	double _mutation_probability;
	
	/* Data Tracking */
	int _do_track;
	vector<Organism*> _initial_organisms;
	vector<Organism*> _current_organisms;
	vector<double> _evaluations;
	
	
	/* PRIVATE METHODS */
	/* Evolutionary Algorithms */
	Manager *alg_1(int run_num);
	Manager *alg_2(int run_num);
	Manager *alg_3(int run_num);
	
	/* Data Tracking Helpers */
	void phylogeny_to_file( std::ofstream& file , Organism* org_ref , int gen );
	void assign_index( Organism* org_ref , int& index );
	void create_lineages( Organism* org_ref , std::vector<Lineage*>& lineages , int generation );
	
	/* Destructor Helpers */
	void delete_organism( Organism* org_ref );
};

#endif
