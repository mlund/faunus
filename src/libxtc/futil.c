/*
 * $Id: /local/xtcio/src/futil.c 21 2006-08-05T20:58:43.457505Z svm  $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "futil.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "fatal.h"
#include "smalloc.h"

bool fexist(const char *fname)
{
  FILE *test;
  
  if (fname == NULL)
    return FALSE;
  test=fopen(fname,"r");
  if (test == NULL) {
    return FALSE;
  } else {
    fclose(test);
    return TRUE;
  }
}

const char *backup_fn(const char *file)
{
  /* Use a reasonably low value for countmax; we might
   * generate 4-5 files in each round, and we dont
   * want to hit directory limits of 1024 or 2048 files.
   */
#define COUNTMAX 128
  int         i,count=1;
  char        *directory,*fn;
  static char buf[256];
  
  for(i=strlen(file)-1; ((i > 0) && (file[i] != '/')); i--)
    ;
  /* Must check whether i > 0, i.e. whether there is a directory
   * in the file name. In that case we overwrite the / sign with
   * a '\0' to end the directory string .
   */
  if (i > 0) {
    directory    = strdup(file);
    directory[i] = '\0';
    fn           = strdup(file+i+1);
  }
  else {
    directory    = strdup(".");
    fn           = strdup(file);
  }
  do {
    sprintf(buf,"%s/#%s.%d#",directory,fn,count);
    count++;
  } while ((count < COUNTMAX) && fexist(buf));
  
  /* Arbitrarily bail out */
  if (count == COUNTMAX) 
    gmx_fatal(FARGS,"Won't make more than %d backups of %s for you",
	      COUNTMAX,fn);
  
  sfree(directory);
  sfree(fn);
  
  return buf;
}
