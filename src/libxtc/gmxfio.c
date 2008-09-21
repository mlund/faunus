/*
 * $Id: /local/xtcio/src/gmxfio.c 27 2008-09-14T05:47:11.743296Z lidb  $
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

#include "gmxfio.h"

#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include "fatal.h"
#include "smalloc.h"
#include "typedefs.h"
#include "macros.h"
#include "futil.h"
#include "xdrf.h"

typedef struct {
  bool bOpen,bRead,bDouble;
  const char *fn;
  XDR  *xdr;
} t_fileio;

static t_fileio *FIO = NULL;
static int  nFIO = 0;

#define fio_check(fio) range_check(fio,0,nFIO)

/*****************************************************************
 *
 *                     EXPORTED SECTION
 *
 *****************************************************************/
int fio_open(const char *fn,
             const char *mode)
{
  t_fileio *fio=NULL;
  int      i,nfio=0;
  const char* bf;
  char     newmode[5];

  if (mode[0]=='r')
    strcpy(newmode,"r");
  else if (mode[0]=='w')
    strcpy(newmode,"w");
  else if (mode[0]=='a')
    strcpy(newmode,"a");
  else
    gmx_fatal(FARGS,"DEATH HORROR in fio_open, mode is '%s'",mode);
  strcat(newmode,"b");
  /* Determine whether we have to make a new one */
  for(i=0; (i<nFIO); i++)
    if (!FIO[i].bOpen) {
      fio  = &(FIO[i]);
      nfio = i;
      break;
    }
  if (i == nFIO) {
    nFIO++;
    srenew(FIO,nFIO);
    fio  = &(FIO[nFIO-1]);
    nfio = nFIO-1;
  }

  fio->bRead  = (newmode[0]=='r');
  fio->bDouble= (sizeof(real) == sizeof(double));
  fio->bOpen  = TRUE;
  fio->fn = strdup(fn);
  fio->xdr = NULL;
  if (newmode[0]=='w') {
    if (fexist(fn)) {
      bf=backup_fn(fn);
      if (rename(fn,bf) == 0) {
        fprintf(stderr,"\nBack Off! I just backed up %s to %s\n",fn,bf);
      } else {
        fprintf(stderr,"Sorry, I couldn't backup %s to %s\n",fn,bf);
      }
    }
  } else {
    /* Check whether file exists */
    if (!fexist(fn))
      gmx_open(fn);
  }
  snew(fio->xdr,1);
  if (!xdropen(fio->xdr,fn,newmode))
    gmx_open(fn);

  return nfio;
}


void fio_close(int fio)
{
  fio_check(fio);
  
  xdrclose(FIO[fio].xdr);
  sfree(FIO[fio].xdr);
  sfree((void *)FIO[fio].fn);
  FIO[fio].bOpen = FALSE;
}

XDR *fio_getxdr(int fio)
{
  fio_check(fio);
  if (FIO[fio].xdr) 
    return FIO[fio].xdr;
  else
    return NULL;
}
