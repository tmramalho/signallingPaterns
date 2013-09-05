/*
 * Signalling Genetics
 * --------------
 * Michael Celentano
 */
 
#include <iostream>
#include <fstream>
#include <istream>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iterator>
# include <set>

/* Boost Random Libraries */
# include <boost/random/uniform_int_distribution.hpp>
# include <boost/random/uniform_real_distribution.hpp>
# include <boost/random/mersenne_twister.hpp>

/* Other Boost Libraries */
# include <boost/lexical_cast.hpp>

#include "Manager.h"
#include "nexist.cpp"
#include "ConstructionMethod.h"
#include "Operations.h"
#include "old/opt/FitnessFunction.h"
#include "old/opt/DNScore.h"
#include "old/opt/TempFF.h"
#include "old/genetic/Evolution.h"
#include "old/genetic/EvAlgorithm.h"
#include "old/genetic/Organism.h"

#include "Reactions/Reaction.h"
#include "Reactions/HillPromReaction.h"
#include "Reactions/HillRepReaction.h"
#include "Reactions/DegReaction.h"

#include "Molecules/Gene.h"
#include "Molecules/Molecule.h"

#include "Operations.h"


using namespace std;

/* This strange-looking pragma is here to disable a warning from Visual C++
 * about truncating long identifiers for debugging symbols. The warning is
 * harmless, but a little disconcernting, so we suppress it. It comes up
 * using STL and other long template expansions.
 */
#if defined(_MSC_VER)
#pragma warning(disable: 4786)
#endif

/*
 * Function macro: main
 * --------------------
 * The purpose of this macro definition is to rename the student
 * main to Main in order to allow a custom main defined in our
 * libraries to configure the application before passing control
 * back to the student program.
 */

#ifdef __APPLE__
#define main Main
#endif

void collect_algorithm_data();

int main () {
	
	
	
	return 0;
}

void collect_algorithm_data() {
	DNScore ff;
	std::string filename;
	std::ofstream file;
	EvAlgorithm alg_type;
	
	for (int i = 0; i < 3; i++) {
		if (i == 0) alg_type = ALG1;
		if (i == 1) alg_type = ALG2;
		if (i == 2)	alg_type = ALG3;
		
		for (int j = 0; j < 10; j++) {
			
			Evolution *ev = new Evolution(&ff);
			ev->colonize(TWO_PROTEIN);
			Manager* best_ref = ev->train(alg_type,i+1);
			
			vector<Lineage*> lineages = ev->create_lineages();
			
			filename = "alg";
			filename += boost::lexical_cast<std::string>(i+1);
			filename += "/";
			filename += "run";
			filename += boost::lexical_cast<std::string>(j+1);
			filename += "/phylogeny.txt";
			file.open(filename.c_str());
			ev->phylogeny_to_file(file);
			file.close();
			
			filename = "alg";
			filename += boost::lexical_cast<std::string>(i+1);
			filename += "/";
			filename += "run";
			filename += boost::lexical_cast<std::string>(j+1);
			filename += "/best_genome.txt";
			file.open(filename.c_str());
			best_ref->genome_to_file(file,"\t");
			file.close();
			
			filename = "alg";
			filename += boost::lexical_cast<std::string>(i+1);
			filename += "/";
			filename += "run";
			filename += boost::lexical_cast<std::string>(j+1);
			filename += "/all_lineages.txt";
			file.open(filename.c_str());
			for ( int k = 0 ; k < lineages.size() ; k++ ) {
				lineages.at(k)->lineage_to_file(file);
			}
			file.close();
			
			filename = "alg";
			filename += boost::lexical_cast<std::string>(i+1);
			filename += "/";
			filename += "run";
			filename += boost::lexical_cast<std::string>(j+1);
			filename += "/evaluations.txt";
			file.open(filename.c_str());
			ev->evaluations_to_file(file);
			file.close();
			
			delete ev;
			for (int k = 0; k<lineages.size(); k++) {
				delete lineages.at(k);
			}
			lineages.clear();
		}
	}	
}


