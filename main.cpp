/*
 * Signalling Genetics
 * --------------
 * Michael Celentano
 */
 
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <vector>

/* Boost Random Libraries */
# include <boost/random/uniform_int_distribution.hpp>
# include <boost/random/uniform_real_distribution.hpp>
# include <boost/random/mersenne_twister.hpp>

#include "Manager.h"
#include "nexist.cpp"
#include "ConstructionMethod.h"
#include "old/opt/FitnessFunction.h"
#include "old/opt/DNScore.h"
#include "old/opt/TempFF.h"
#include "old/genetic/Evolution.h"

#include "Reactions/Reaction.h"
#include "Reactions/PromReaction.h"
#include "Reactions/DegReaction.h"


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

/*Main program*/

int main() 
{

	DNScore ff;
	Evolution ev(&ff);
	ev.colonize(TWO_PROTEIN);
	
	Manager *best_ref = ev.train();
	
	std::ofstream file;
	file.open("best_genome.txt");
	best_ref->genome_to_file(file,"\t");
	file.close();
	file.open("best_integration.txt");
	best_ref->initialize();
	unsigned seed = std::time(0);
	boost::random::mt19937 generator(seed);
	best_ref->integrate(10000,generator,file);
	file.close();
	
	/*
	for ( int i = 0 ; i < sc_ref->_neighbors.size() ; i++ ) {
		for ( int j = 0 ; j < sc_ref->_neighbors.at(i)->size() ; j++ ) {
			std::cout << sc_ref->_neighbors.at(i)->at(j) << " :::: ";
		}
		std::cout << std::endl;
	}
	*/
	
	
	
	
	/*
	SettingsCont* sc_ref = SettingsCont::getInstance();
	sc_ref->set_neighbors("1,|0,2,|1,3,|2,4,|3,5,|4,6,|5,7,|6,8,|7,9,|8,10,|9,11,|10,12,|11,13,|12,14,|13,15,|14,16,|15,17,|16,|");
	Manager manager(COLLIER_DELTA_NOTCH);
	
	std::ofstream file;
	file.open("collier_delta_notch_18_stc_point05_trial.txt");
	unsigned seed = std::time(0);
	boost::random::mt19937 generator(seed);
	manager.integrate(30000, generator, file);
	file.close();
	file.open("collier_delta_notch_genome_stc.txt");
	manager.genome_to_file(file,"\t");
	file.close();
	
	/*
	TempFF ff;
	Evolution evolution(&ff);
	evolution.colonize("one_protein");
	Manager* best = evolution.train();
	
	std::ofstream file;
	file.open("genome.txt");
	best->genome_to_file(file,"\t");
	file.close();
	file.open("history.txt");
	best->initialize();
	unsigned seed = std::time(0);
	boost::random::mt19937 generator(seed);
	best->integrate(1000,generator,file);
	file.close();
	 */
	
	return 0;
	
}