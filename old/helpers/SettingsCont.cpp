/*
 * SettingsCont.cpp
 *
 *  Created on: Jan 12, 2012
 *      Author: tiago
 */

#include "SettingsCont.h"

SettingsCont *SettingsCont::_sc_ref = NULL;

SettingsCont::SettingsCont() {
	_num_threads = 2;
	_verbose = 0;
	_na = 20;
	_ng = 50;
	_stochasticity = .01;
	
	_neighbors.push_back(new std::vector<int>(1,1));
	std::vector<int>* addition = new std::vector<int>;
	addition->push_back(0);
	addition->push_back(2);
	_neighbors.push_back(addition);
	addition = new std::vector<int>(1,1);
	_neighbors.push_back(addition);
	
	_mode = RK4_STC_TI;
	_dt = .001;
}
/*
void SettingsCont::setParameters(int argc, char *argv[]) {
	int c;
	while ((c = getopt (argc, argv, "v:t:d:c:p:s:n:a:u:g:b:")) != -1) {
		switch (c) {
			case 'g':
				ng = atoi(optarg);
				break;
			case 'u':
				na = atoi(optarg);
				break;
			case 'a':
				acc = atof(optarg);
				break;
			case 't':
				numThreads = atoi(optarg);
				break;
			case 'v':
				verbose = atoi(optarg);
				break;
			case 'd':
				dim = atoi(optarg);
				break;
			case 'b':
				bound = atof(optarg);
				break;
			case 'c':
				coop = atof(optarg);
				break;
			case 'p':
				p = atof(optarg);
				break;
			case 's':
				sigma0 = atof(optarg);
				break;
			case 'n':
				poiss0 = atof(optarg);
				break;
			case '?':
				if (optopt == 't')
					fprintf (stderr, "Option -%c requires an argument.\n", optopt);
				else if (isprint (optopt))
					fprintf (stderr, "Unknown option `-%c'.\n", optopt);
				else
					fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
			default:
				abort ();
		}
	}
	std::cout << "Parameters have been set:" << std::endl;
	std::cout << std::setw(20) << "-t numThreads: "    << numThreads << std::endl;
	std::cout << std::setw(20) << "-v verbose: "       << verbose << std::endl;
	std::cout << std::setw(20) << "-d dimension: "     << dim << std::endl;
	std::cout << std::setw(20) << "-b minimum conc.: " << bound << std::endl;
	std::cout << std::setw(20) << "-c cooperativity: " << coop << std::endl;
	std::cout << std::setw(20) << "-p p: "             << p << std::endl;
	std::cout << std::setw(20) << "-s \\sigma0: "      << sigma0 << std::endl;
	std::cout << std::setw(20) << "-n N0: "            << poiss0 << std::endl;
	std::cout << std::setw(20) << "-g numGenerations: "<< ng << std::endl;
	std::cout << std::setw(20) << "-u numAgents: "     << na << std::endl;
	std::cout << std::setw(20) << "-a GSS accuracy: "  << acc << std::endl;
	std::cout << std::endl;
}
*/
SettingsCont::~SettingsCont() {
	_neighbors.erase(_neighbors.begin(),_neighbors.end());
}

SettingsCont *SettingsCont::getInstance() {
	if(_sc_ref == NULL) {
		_sc_ref = new SettingsCont();
	}

	return _sc_ref;
}

void SettingsCont::set_neighbors ( std::string neighbors_code ) {
	
	_neighbors.erase(_neighbors.begin(),_neighbors.end());
	
	
	for ( std::string::iterator it = neighbors_code.begin() ;
		 it != neighbors_code.end() ;
		 it++ ) {
	
		std::vector<int>* next = new std::vector<int>;
		
		/* While we are still reading this in this vector */
		while ( *it != '|' ) {
			
			std::string num = "";
			/* While we are still reading this number */
			while ( *it != ',' ) {
				num += *it;
				it++;
			}
			next->push_back(atoi(num.c_str()));
			it++;
			
		}
		_neighbors.push_back(next);
		
	}
	
}



