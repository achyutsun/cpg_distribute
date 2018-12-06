#include <stdio.h>

/********************* CONTROLS ********************/
#define CPGSCORE  17   /* Changed from 28 960214 */
/* so that can compare with old cpg - this
			  had CPGSCORE 27, but allowed score to reach
			  0 once without reporting */

void process_fasta(FILE *f);
