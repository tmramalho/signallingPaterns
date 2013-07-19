/*
 *  Gene.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/17/13.
 *
 */

#include "Gene.h"

Gene::Gene( int iProduct , int iBoundPromoter , int iRoot ) {
	this->iProduct = iProduct;
	this->iBoundPromoter = iBoundPromoter;
	this->iRoot = iRoot;
}

/* Destructor: ~Gene()
 * -------------------------------------------------------------------------- 
 * No heap allocated memory owned by the gene class.
 */
Gene::~Gene() {}