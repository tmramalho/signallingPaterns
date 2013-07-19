/*
 *  Protein.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/17/13.
 *
 */

#include "Protein.h"

/* Constructor: Protein()
 * -------------------------------------------------------------------------- 
 */
Protein::Protein( int iSelf , int iRootZero , int iRootOne ) {
	this->iRootZero = iRootZero;
	this->iRootOne = iRootOne;
	
	i = iSelf;
	
	/* reactions will automatically initiate to empty vector */
}

/* Destructor: ~Protein()
 * -------------------------------------------------------------------------- 
 * No heap allocated memory owned by the protein class.
 */
Protein::~Protein() {}
