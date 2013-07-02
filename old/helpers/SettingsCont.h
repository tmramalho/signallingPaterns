/*
 * SettingsCont.h
 *
 *  Created on: Jan 12, 2012
 *      Author: tiago
 *
 *  Singleton design pattern
 */

#ifndef SETTINGSCONT_H_
#define SETTINGSCONT_H_

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <unistd.h>

#define SIG(x,b) ((1.0-b)/(1.0+exp(-x))+b)

class SettingsCont {

private:
	static SettingsCont *scPointer;
	SettingsCont();

public:
	int numThreads, dim, na, ng;
	int verbose;
	double bound, coop, p, sigma0, poiss0, acc;
	static SettingsCont *getInstance();
	void setParameters(int argc, char *argv[]);
	~SettingsCont();
};

#endif /* SETTINGSCONT_H_ */
