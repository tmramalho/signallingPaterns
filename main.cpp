/*
 * Signalling Genetics
 * --------------
 * Michael Celentano
 */
 
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "ODEManager.h"


/*Main program*/

int main() 
{
	ODEManager manager;
	manager.run("rk1_det_ti");
	
	return 0;
}