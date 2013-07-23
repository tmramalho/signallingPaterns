/*
 *  Gene.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/17/13.
 *
 */

#include "Gene.h"

Gene::Gene( int iSelf, int iProduct , int iBoundPromoter , int iRoot ) {
	this->iProduct = iProduct;
	this->iBoundPromoter = iBoundPromoter;
	this->iRoot = iRoot;
	i = iSelf;
}

/* Destructor: ~Gene()
 * -------------------------------------------------------------------------- 
 * No heap allocated memory owned by the gene class.
 */
Gene::~Gene() {}