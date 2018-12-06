#include "options.h"
#include "cpg.h"
#include <stdio.h>

/*
==============================================================================
Copyright (c) 2004  Gos Micklem + Tim Cutts

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

Contact:  g.micklem@gen.cam.ac.uk
==============================================================================
*/

int main (int argc, char *argv[]) { 
  
    FILE *fil;

    int last_arg;
    
    last_arg = options(argc, argv);

    if (argc == last_arg)
        fil = stdin;
    else if (!(fil = fopen ( argv[last_arg], "r" ))) {
	perror("Could not open file");
	return 1;
    }
  
    process_fasta(fil);

    fclose(fil);

    return 0;
}
