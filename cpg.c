/*  $Id: cpg.c,v 1.2 2004/03/04 20:17:01 tim Exp $  */
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

  cpg.c: CpG island finder. Larsen et al Genomics 13 1992 p1095-1107
  "usually defined as >200bp with %GC > 50% and obs/exp CpG >
  0.6". Here use running sum rather than window to score: if not CpG
  at postion i, then decrement runSum counter, but if CpG then runSum
  += CPGSCORE.     Spans > threshold are searched
  for recursively.
  
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "cpg.h"
#include "options.h"

void getstats (int start, int end, char *seq, char *seqname,
		int *ncpg, int *ngpc, int *ngc) {

    int i;
    *ncpg = *ngpc = *ngc = 0;
    
	/*

	Code used to look like this.  Newer version looks more
	complicated but has far fewer comparisons during execution.

    for (i = start; i < end; ++i)  {
	if (end-1-i && seq[i]=='C' && seq[i+1]=='G') ++*ncpg; 
	if (end-1-i && seq[i]=='G' && seq[i+1]=='C') ++*ngpc; 
	if (seq[i]=='C' || seq[i]=='G') ++*ngc; 
	*/
    
    end--;

    for (i = start; i < end; ++i) {
	if (seq[i]=='C') {
	    ++*ngc;
	    if (seq[i+1]=='G')
		++*ncpg;
	} else if (seq[i]=='G') {
	    ++*ngc;
	    if (seq[i+1]=='C')
		++*ngpc;
	}
    }

    if (seq[end]=='C' || seq[end]=='G')
	++*ngc;

}

void print_output(char *seqname, int start, int stop, int score,
		  int ncpg, double gc_pct, double oe) {

    if (stop-start < length_threshold)
	return;

    if (gc_pct < gc_threshold)
	return;

    if (oe < observed_expected_threshold)
	return;

    if (ace_output) {
	printf("Feature Predicted_CpG_island %d %d %.2f ",
	       start, stop, oe);
	printf("\"Predicted CpG island: %.1f %%GC,   o/e=%.2f,   #CpGs= %d\"\n",
	       gc_pct, oe, ncpg);
    } else {
	printf("%s\t%d\t%d\t%d\tCpG:%d\t%.1f\t",
	       seqname, start, stop, score,
	       ncpg, gc_pct);
	
	if (oe >= 0.0)
	    printf("%.2f\n", oe);
	else
	    puts("-");
    }
}


void findspans (int start, int end, char *seq, char *seqname){
    int i ;
    int sc = 0;
    int lsc = 0; 
    int imn = -1;  /* else sequences which start with pattern reported badly */
    int imx = 0;
    int mx = 0;
    
    int ncpg, ngpc, ngc;  
    
    i = start;
    while (i < end)  {
	lsc = sc;
	if (seq[i] == 'C' && seq[i+1] == 'G' && end-1-i)
	    sc += CPGSCORE;
	else {
	    if (sc > 0)
		--sc;
	}
	
	/*      printf("%d \t %d \t%d \t %d \t%d \t%d\n", i, sc, lsc, imn, imx, mx); */
	if (sc == 0 && lsc) {
	    /* imn+1 is the start of the match. 
	       imx+1 is the end of the match, given pattern-length=2.
	       fetchdb using reported co-ordinates will return a sequence
	       which both starts and ends with a complete pattern.
	       Further +1 adjusts so that start of sequence = 1 
	    */
	    
	    getstats(imn+1, imx+2, seq, seqname, &ncpg, &ngpc, &ngc);

	    print_output(seqname, imn+2, imx+2, mx, ncpg,
			 ngc*100.0/(imx+1-imn),
			 ngpc == 0 ? -1 : 1.0*ncpg/ngpc); 

	    /* 	  printf("%s \t %d\t %d\t %d \n", seqname, imn+2, imx+2, mx); 
	     */
	    /* Recursive call searches from one position after the end of the 
	       last match to the current position */
	    findspans(imx+2, i, seq, seqname);
	    sc = lsc = imn = imx =  mx = 0;
	}
	imx = sc > mx ? i : imx;
	mx = sc > mx ? sc : mx;
	imn = sc == 0 ? i : imn;
	++i;
    }

    if (sc != 0)  {
	/*      printf("%d \t %d \t%d \t %d \t%d \t%d\n", i, sc, lsc, imn, imx, mx);  */
	
	getstats (imn+1, imx+2, seq, seqname, &ncpg, &ngpc, &ngc);
	
	print_output(seqname, imn+2, imx+2, mx, ncpg,
		     ngc*100.0/(imx+1-imn),
		     ngpc == 0 ? -1 : 1.0*ncpg/ngpc); 

	/*      printf("%s \t %d\t %d\t %d \n", seqname, imn+2, imx+2, mx); 
	 */
	findspans(imx+2, end, seq, seqname);
    }
}

void process_sequence(int len, char *s, char *title) {
    char *p;
    
    /* Force to upper case */
    for (p = s; *p != '\0'; p++) *p &= ~0x20;
    
    /* Truncate title to sequence ID */
    for (p = title; !isspace(*p); p++) ;
    
    *p = '\0';
    
    findspans(0, len, s, title);
}

static int sz=8192;

void process_fasta(FILE *f)
{
    char *s, *p;
    char buf[1024];
    char title[1024];
    int n, len;
    
    s = (char *)malloc(sz);
    if (!s) {
	fprintf(stderr, "Out of memory\n");
	exit(1);
    }

    *s = '\0'; len = 0;

    while (!feof(f)) {
	p = fgets(buf, sizeof(buf), f);
	if (p) {
	    if (*p == '>') {
		/* We have found a comment line, so we can process */
		strcpy(title, ++p);
		if (*s != '\0') {
		    s[len] = '\0';

		    process_sequence(len, s, title);

		    len = 0; *s= '\0';
		}
	    } else {
		n = strlen(p)-1;
		if ((len+n+3)>sz) {
		    sz <<= 1;
		    s = (char *)realloc(s, sz);
		    if (!s) {
			fprintf(stderr, "Out of memory\n");
			exit(1);
		    }
		}
		memcpy(&s[len], buf, n);
		len += n;
	    }
	}
    }

    if (*s != '\0') {
	s[len] = '\0';

	process_sequence(len, s, title);

	len = 0; *s= '\0';
    }
}
