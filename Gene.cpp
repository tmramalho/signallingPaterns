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

Gene::~Gene() {}