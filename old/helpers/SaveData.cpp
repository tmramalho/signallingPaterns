/*
 * SaveData.cpp
 *
 *  Created on: Jan 16, 2012
 *      Author: tiago
 */

#include "SaveData.h"

SaveData::SaveData() {
	// TODO Auto-generated constructor stub

}

SaveData::~SaveData() {
	// TODO Auto-generated destructor stub
}

void SaveData::saveMathematica(dvec *data) {
	Profile *m = new Profile(data);
	SettingsCont *s = SettingsCont::getInstance();
	std::fstream filestr("mPoints.txt", std::fstream::app | std::fstream::out);

	filestr << "{";
	for(int i=0; i<100; i++) {
		double x =(double)i/(double)(100-1);
		filestr << "{ " << x << " , " << SIG(m->eval(x),s->bound) << " }";
		if(i!=99) filestr << ",";
		filestr << std::endl;
	}
	filestr << "}";
	filestr.close();
}

void SaveData::saveData(dvec *data, double sc) {

	Profile *m = new Profile(data);
	SettingsCont *s = SettingsCont::getInstance();
	char buffer[256];
	sprintf(buffer,"out%ld_%d.log", time(NULL), getpid());
	std::fstream filestr(buffer, std::fstream::app | std::fstream::out);

	filestr << "Final score:" << sc << std::endl;
	filestr << "Final values:" << *data << std::endl << std::endl;
	filestr << "Parameters: " << std::endl;
	filestr << std::setw(20) << "-t numThreads: "    << s->numThreads << std::endl;
	filestr << std::setw(20) << "-v verbose: "       << s->verbose << std::endl;
	filestr << std::setw(20) << "-d dimension: "     << s->dim << std::endl;
	filestr << std::setw(20) << "-b minimum conc.: " << s->bound << std::endl;
	filestr << std::setw(20) << "-c cooperativity: " << s->coop << std::endl;
	filestr << std::setw(20) << "-p p: "             << s->p << std::endl;
	filestr << std::setw(20) << "-s \\sigma0: "      << s->sigma0 << std::endl;
	filestr << std::setw(20) << "-n N0: "            << s->poiss0 << std::endl;
	filestr << std::setw(20) << "-g numGenerations: "<< s->ng << std::endl;
	filestr << std::setw(20) << "-u numAgents: "     << s->na << std::endl;
	filestr << std::setw(20) << "-a GSS accuracy: "  << s->acc << std::endl;
	filestr << std::endl;

	filestr.close();

	//csv with the spline
	sprintf(buffer,"out%ld_%d.csv", time(NULL), getpid());
	std::fstream csvstr(buffer, std::fstream::app | std::fstream::out);

	int imax = 500;
	for(int i=0; i<imax; i++) {
		double x =(double)i/(double)(imax-1);
		csvstr << x << " , " << m->eval(x) << std::endl;
	}
	csvstr.close();
}
