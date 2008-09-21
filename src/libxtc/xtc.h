/*
 * $Id: /local/xtcio/src/xtc.h 30 2008-09-14T06:02:45.605229Z lidb  $
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

#ifndef XTC_H
#define XTC_H

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <rpc/types.h>
#include <rpc/xdr.h>

#if defined(XTC_ENABLE_FLOAT)
typedef float xtc_real;
#else
typedef double xtc_real;
#endif

typedef xtc_real xtc_rvec[3];
typedef xtc_real xtc_matrix[3][3];

typedef struct Xtc Xtc;
#define XTC_MAGIC 1995

#ifdef __cplusplus
extern "C" {
#endif

Xtc* xtc_create(const char* filename, const char* mode);
Xtc* xtc_stdio_create(FILE* file, const char* mode);
void xtc_destory(Xtc* self);

bool xtc_read_first(Xtc* self,
                   int *natoms,
                   int *step,
                   xtc_real *time,
                   xtc_matrix box,
                   xtc_rvec **x,
                   xtc_real *prec);
/* Open xtc file, read xtc file first time, allocate memory for x */

bool xtc_read_next(Xtc* self,
                  int natoms,
                  int *step,
                  xtc_real *time,
                  xtc_matrix box,
                  xtc_rvec *x,
                  xtc_real *prec);
/* Read subsequent frames */

bool xtc_write(Xtc* self,
              int natoms,
              int step,
              xtc_real time,
              xtc_matrix box,
              xtc_rvec *x,
              xtc_real prec);
/* Write a frame to xtc file */


bool
xtc_seek_time(Xtc* self, 
              double time, 
              int natoms);

bool
xtc_seek_frame(Xtc* self,
               int frame, 
               int natoms);

double
xtc_get_last_frame_time(Xtc* self,
                        int natoms, 
                        bool* bOK);

int 
xtc_get_last_frame_number(Xtc* self, 
                          int natoms, 
                          bool* bOK);

#ifdef __cplusplus
}
#endif

#endif
