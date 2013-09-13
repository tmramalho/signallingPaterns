/*
 * Signalling Genetics
 * --------------
 * Michael Celentano
 */
 
#include "cs106/scanner.h"

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

#include "ManagerDebugger.h"


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
void reaction_deletion_debugging();
void collier_delta_notch_run();
void test_genome_in_file();

void create_documentation( ifstream& input );

int main () {
	
	ifstream input;
	input.open("DocumentationCode.txt");
	create_documentation(input);
	
	return 0;
}

/*
void collect_algorithm_data() {
	DNScore ff;
	std::string filename;
	std::ofstream file;
	EvAlgorithm alg_type;
	
	for (int i = 2; i < 3; i++) {
		if (i == 0) alg_type = ALG1;
		if (i == 1) alg_type = ALG2;
		if (i == 2)	alg_type = ALG3;
		
		for (int j = 5; j < 10; j++) {
			
			Evolution *ev = new Evolution(&ff);
			ev->colonize(TWO_PROTEIN);
			Manager* best_ref = ev->train(alg_type,j+1);
			
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
*/
void reaction_deletion_debugging() {
	ManagerDebugger md;
	md.reaction_removal_debugging();
}

void collier_delta_notch_run() {
	Manager *m = new Manager(COLLIER_DELTA_NOTCH);
	m->initialize();
	ofstream file;
	file.open("genome.txt");
	m->genome_to_file(file,"\t");
	file.close();
	file.open("integration.txt");
	unsigned seed = std::time(0);
	boost::random::mt19937 generator(seed);
	m->integrate(10000,generator,file);
	file.close();
}

void test_genome_in_file() {
	
	unsigned seed = std::time(0);
	boost::random::mt19937 generator(seed);
	
	ifstream file;
	file.open("alg2/run2/best_genome.txt");
	Manager *m = new Manager(file);
	file.close();
	
	ofstream output;
	output.open("alg2/run2/best_genome_test2.txt");
	m->initialize();
	m->integrate(50000, generator, output);
	output.close();
	
}

void create_documentation( ifstream& input ) {
	
	int num_indents = 0;
	int counter = 0;
	int _line_size = 75;
	
	ofstream output;
	output.open("Documentation.txt");
	
	Scanner scan;
	scan.setInput(input);
	scan.setSpaceOption(Scanner::PreserveSpaces);
	scan.setNumberOption(Scanner::ScanNumbersAsIntegers);
	scan.setStringOption(Scanner::ScanQuotesAsPunctuation);
	scan.setBracketOption(Scanner::ScanBracketsAsPunctuation);
	
	string next_token;
	
	while (scan.hasMoreTokens()) {
		next_token = scan.nextToken();
		if (next_token == "LINEBREAK") {
			output << "\n";
			for (int i = 0; i < num_indents; i++) {
				output << "\t";
			}
			counter = 7 * num_indents;
		}
		else if ( next_token == "NEWPARAGRAPH" ) {
			output << "\n" << "\n";
			for (int i = 0; i < num_indents; i++) {
				output << "\t";
			}
			counter = 7 * num_indents;
		}
		else if (next_token == "SETINDENT") {
			scan.nextToken();
			num_indents = atoi(scan.nextToken().c_str());
		}
		else if (next_token == "INCREMENTINDENT") {
			scan.nextToken();
			string next = scan.nextToken();
			if ( next == "-" ) {
				num_indents -= atoi(scan.nextToken().c_str());
			}
			else {
				num_indents += atoi(next.c_str());
			}

		}
		else if ( next_token == "CROSSLINE") {
			output << "\n";
			for (int i = 0; i < _line_size; i++) {
				output << '-';
			}
			output << "\n";
			for (int i = 0; i < num_indents; i++) {
				output << "\t";
			}
			counter = 7 * num_indents;
		}
		else if (next_token == "\n" || next_token == "\t");
		else {
			counter += next_token.size();
			string next = scan.nextToken();
			if ( next_token != " " && next != " " && next != "\t" && next != "\n") {
				while (next != " " && next != "\t" && next != "\n") {
					next_token += next;
					counter += next.size();
					next = scan.nextToken();
				}
			}
			scan.saveToken(next);
			if ( counter >= _line_size ) {
				output << "\n";
				for (int i = 0; i < num_indents ; i++) {
					output << "\t";
				}
				counter = 7 * num_indents;
				if ( next_token != " " ) {
					output << next_token;
					counter += next_token.size();
				}
			}
			else {
				output << next_token;
			}
		}
	}
	
	output.close();
	
}
