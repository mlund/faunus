/*
 * $Id: /local/xtcio/src/xdrf.h 39 2008-09-14T06:31:00.970588Z lidb  $
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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _xdrf_h
#define _xdrf_h

#include <stdio.h>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include <limits.h>
#include "xtc.h"

/* Read or write reduced precision *float* coordinates */
int 
xdr3dfcoord(XDR *xdrs, 
            float *fp, 
            int *size, 
            float *precision);

int xdropen(XDR *xdrs, const char *filename, const char *type);
int xdrclose(XDR *xdrs);
float _xtc_estimate_dt(int fp, int natoms, bool * bOK);
float _xtc_get_next_frame_time(int fp,int natoms, bool * bOK);

#ifdef HAVE_FSEEKO
#define gmx_off_t off_t
#else
#define gmx_off_t int
#endif

gmx_off_t _xtc_get_next_frame_start(int fp, int natoms);
int _xtc_get_next_frame_number(int fp,int natoms);
float _xtc_get_current_frame_time(int fp,int natoms, bool * bOK);
int _xtc_get_current_frame_number(int fp,int natoms, bool * bOK);
extern FILE* xdrfiles[];

#endif







