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
	//std::cout << "Reaction Default Constructor\n";
	_sc_ref = SettingsCont::getInstance();
}

Reaction::~Reaction() {
	//std::cout << "Reaction Destructor " << _type << "\n";
}

bool Reaction::operator==( const Reaction& rhs ) const {
	
	if ( _type != rhs._type ) return false;
	else {
		for( int part_num = 0 ; part_num < get_num_part() ; part_num++ ) {
			if ( get_i_part(part_num) != rhs.get_i_part(part_num) ) return false;
		}
	}
	return true;
	
}

bool Reaction::operator!=( const Reaction& rhs ) const {
	return !(*this == rhs);
}

/* Public Method: get_num_part()
 * --------------------------------------------------------------------------
 * Returns the number of participants in this reaction. 
 */
int Reaction::get_num_part() const {
	return _num_part;
}

/* Public Method: get_type()
 * --------------------------------------------------------------------------
 * Returns the type of the reaction.
 */
ReactionType Reaction::get_type() {
	return _type;
}