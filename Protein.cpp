/*
 *  Protein.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/17/13.
 *
 */

#include "Protein.h"

Protein::Protein( int iSelf , int iRootZero , int iRootOne ) {
	this->iRootZero = iRootZero;
	this->iRootOne = iRootOne;
	
	i = iSelf;
	
	/* reactions will automatically initiate to empty vector */
}

Protein::~Protein() {}
