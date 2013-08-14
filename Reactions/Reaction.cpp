/*
 *  Reaction.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/18/13.
 *
 */

#include "Reaction.h"

/* Constructor: Reaction()
 * --------------------------------------------------------------------------
 */
Reaction::Reaction() {
	_sc_ref = SettingsCont::getInstance();
}

/* Public Method: get_num_part()
 * --------------------------------------------------------------------------
 * Returns the number of participants in this reaction. 
 */
int Reaction::get_num_part() {
	return _num_part;
}

/* Public Method: get_type()
 * --------------------------------------------------------------------------
 * Returns the type of the reaction.
 */
ReactionType Reaction::get_type() {
	return _type;
}