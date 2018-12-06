#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <errno.h>
#include "options.h"

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

int length_threshold = DEF_LENGTH_THRESHOLD;
double gc_threshold = DEF_GC_THRESHOLD;
double observed_expected_threshold = DEF_OBSEXP_THRESHOLD;
int ace_output = 0;

void usage (char *prog) { 
    fprintf(stderr, "\nUsage: %s [options] [file]\n\n", prog);
    fprintf(stderr, "\tOptions:\n\n");
    fprintf(stderr, "\t\t-h\tThis message\n\n");
    fprintf(stderr, "\t\t-a\t.ace format output\n\n");
    fprintf(stderr, "\t\t-l nn\tlength threshold\n");
    fprintf(stderr, "\t\t-g nn\t%%GC threshold\n");
    fprintf(stderr, "\t\t-o nn\tobserved/expected threshold\n\n");
    fprintf(stderr, "\tStringent CpG island: -l 1000 -g 50 -o 0.6\n\n");
}

int options (int argc, char *argv[]) {

    int c;
    int errflg = 0;
    char *endptr;

    while (!errflg && (c = getopt(argc, argv, "ahl:g:o:")) != -1) {
	switch (c) {
	case 'h': usage(argv[0]); exit(0); break;
	case 'a': ace_output++; break;
	case 'l': {
	    length_threshold = strtol(optarg, &endptr, 0);
	    if (length_threshold == 0 && *endptr != '\0') {
		fprintf(stderr,
			"%s: Invalid characters found in length threshold: %s\n\n",
			argv[0],
			endptr);
		errflg++;
	    } else if (length_threshold < 0 || length_threshold == LONG_MAX) {
		fprintf(stderr, "%s: Invalid length threshold\n\n", argv[0]);
		errflg++;
	    }
	    break;
	}
	case 'g': {
	    gc_threshold = strtod(optarg, &endptr);
	    if (errno == ERANGE || gc_threshold < 0.0 ||
		gc_threshold > 100.0 ) {
		fprintf(stderr, "%s: Invalid %%GC threshold\n\n", argv[0]);
		errflg++;
	    }
	    break;
	}
	case 'o': {
	    observed_expected_threshold = strtod(optarg, &endptr);
	    if (errno == ERANGE || observed_expected_threshold < 0.0) {
		fprintf(stderr,
			"%s: Invalid observed/expected threshold\n\n",
			argv[0]);
		errflg++;
	    }
	    break;
	}
	default: {
	    errflg++;
	}
	}
    }

    if (errflg)
	exit(errflg);

    return optind;
}
