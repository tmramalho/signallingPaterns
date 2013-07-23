/*
 * Signalling Genetics
 * --------------
 * Michael Celentano
 */
 
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "Manager.h"

using namespace std;

/* This strange-looking pragma is here to disable a warning from Visual C++
 * about truncating long identifiers for debugging symbols. The warning is
 * harmless, but a little disconcernting, so we suppress it. It comes up
 * using STL and other long template expansions.
 */
#if defined(_MSC_VER)
#pragma warning(disable: 4786)
#endif

/*
 * Function macro: main
 * --------------------
 * The purpose of this macro definition is to rename the student
 * main to Main in order to allow a custom main defined in our
 * libraries to configure the application before passing control
 * back to the student program.
 */

#ifdef __APPLE__
#define main Main
#endif

/*Main program*/

int main() 
{
	Manager manager;
	
	manager.integrate(1000);
	
	return 0;
}