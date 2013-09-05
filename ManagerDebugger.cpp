/*
 *  ManagerDebugger.cpp
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 9/5/13.
 *
 */

#include "ManagerDebugger.h"

void ManagerDebugger::reaction_removal_debugging() {
	
	Manager m(TWO_PROTEIN);
	
	/* < a b A B > */
	/* < 0 1 2 3 > */
	
	m.add_gene(0.0,1.0);
	/* < a b c A B C > */
	/* < 0 1 2 3 4 5 > */
	
	m.add_CombReaction(0,1,1.0,1.0);
	/* < a b c A B C A:B > */
	/* < 0 1 2 3 4 5  6  > */
	
	m.add_CombReaction(3,2,1.0,1.0);
	/* < a b c A B C A:B (A:B):C > */
	/* < 0 1 2 3 4 5  6     7    > */
	
	m.add_CombReaction(0,4,1.0,1.0);
	/* < a b c A B C A:B (A:B):C A:((A:B):C) > */
	/* < 0 1 2 3 4 5  6     7        8       > */
	
	/* Reactions vector: 
	 * < promA degA promB degB promC degC comb0 comb1 comb2 > */
	/* <   0     1    2    3     4     5    6     7     8   > */
	
	std::ofstream file;
	file.open("removal_debugging_1.txt");
	m.genome_to_file(file,"\t");
	file.close();
	
	m.remove_reaction ( 6 ) ;
	/* < a b c A B C > */
	/* < 0 1 2 3 4 5 > */
	/* Reactions vector: 
	 * < promA degA promB degB promC degC > */
	/* <   0     1    2    3     4     5  > */
	
	
	file.open("removal_debugging_2.txt");
	m.genome_to_file(file,"\t");
	file.close();
	
	m.add_CombReaction(0,1,1.0,1.0);
	/* < a b c A B C A:B > */
	/* < 0 1 2 3 4 5  6  > */
	
	m.add_CombReaction(3,2,1.0,1.0);
	/* < a b c A B C A:B (A:B):C > */
	/* < 0 1 2 3 4 5  6     7    > */
	
	m.add_CombReaction(0,4,1.0,1.0);
	/* < a b c A B C A:B (A:B):C A:((A:B):C) > */
	/* < 0 1 2 3 4 5  6     7        8       > */
	
	/* Reactions vector: 
	 * < promA degA promB degB promC degC comb0 comb1 comb2 > */
	/* <   0     1    2    3     4     5    6     7     8   > */
	
	file.open("removal_debugging_3.txt");
	m.genome_to_file(file,"\t");
	file.close();
	
	m.remove_reaction(7);
	
	/* < a b c A B C A:B > */
	/* < 0 1 2 3 4 5  6  > */
	/* Reactions vector: 
	 * < promA degA promB degB promC degC comb0 > */
	/* <   0     1    2    3     4     5    6   > */
	
	file.open("removal_debugging_4.txt");
	m.genome_to_file(file,"\t");
	file.close();
	
	m.remove_reaction(2);
	
	/* < a b c A B C A:B > */
	/* < 0 1 2 3 4 > */
	/* Reactions vector: 
	 * < promA degA degB promC degC comb0 > */
	/* <   0     1    2    3    4     5   > */
	
	file.open("removal_debugging_5.txt");
	m.genome_to_file(file,"\t");
	file.close();
	
	m.add_PromBindingReac(1, 3, 1.0, 1.0, 1.0);
	
	/* < a b c b:(A:B) A B C A:B > */
	/* < 0 1 2	 3     4 5 6  7  > */
	/* Reactions vector: 
	 * < promA degA degB promC degC comb0 PromBinding NewPromotion> */
	/* <   0     1    2    3     4     5      6              7    > */
	
	file.open("removal_debugging_6.txt");
	m.genome_to_file(file,"\t");
	file.close();
	
	m.remove_gene(1);
	
	/* < a c A C > */
	/* < 0 1 2 3 > */
	/* Reactions vector: 
	 * < promA degA promC degC > */
	/* <   0     1    2     3  > */
	
	file.open("removal_debugging_7.txt");
	m.genome_to_file(file,"\t");
	file.close();
	
	
}