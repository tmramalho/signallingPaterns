/*
 *  ReactionType.h
 *  Simulating Genetics
 *
 *  Created by Michael Celentano 2 on 7/18/13.
 *
 */

/* ReactionType allows us to classify reactions easily.
 * 
 */

#ifndef REACTIONTYPE_H
#define REACTIONTYPE_H

enum ReactionType { PROMOTION ,
					DEGRADATION , 
					COMBINATION , 
					LATERAL_PROMOTION , 
					PROMOTER_BINDING ,
					HILL_PROMOTION ,
					HILL_REPRESSION };

#endif