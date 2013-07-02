/*
 * SaveData.h
 *
 *  Created on: Jan 16, 2012
 *      Author: tiago
 */

#ifndef SAVEDATA_H_
#define SAVEDATA_H_

#include <fstream>

#include "../opt/FitnessFunction.h"
#include "../numeric/Profile.h"
#include "../helpers/SettingsCont.h"

class SaveData {
public:
	SaveData();
	virtual ~SaveData();

	static void saveMathematica(dvec *data);
	static void saveData(dvec *data, double sc);
};

#endif /* SAVEDATA_H_ */
