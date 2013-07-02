/*
 * SettingsCont.cpp
 *
 *  Created on: Jan 12, 2012
 *      Author: tiago
 */

#include "SettingsCont.h"

SettingsCont *SettingsCont::scPointer = NULL;

SettingsCont::SettingsCont() {
	numThreads = 2;
	verbose = 0;
	dim = 10;
	bound = 0.01;
	coop = 4.0;
	p = 2;
	sigma0 = 0.01;
	poiss0 = 0.01;
	acc = 1e-2;
	na = 20;
	ng = 50;
}

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

SettingsCont::~SettingsCont() {
}

SettingsCont *SettingsCont::getInstance() {
	if(scPointer == NULL) {
		scPointer = new SettingsCont();
	}

	return scPointer;
}
