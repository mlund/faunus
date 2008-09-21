/*
 * $Id: /local/xtcio/src/fatal.h 21 2006-08-05T20:58:43.457505Z svm  $
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

#ifndef FATAL_H
#define FATAL_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

void _where(const char *file, int line);
#define where() _where(__FILE__, __LINE__)
/* Prints filename and line to stdlog and only on amba memvail */

void gmx_fatal(int fatal_errno,const char *file,int line,const char *fmt,...);
#define FARGS 0,__FILE__,__LINE__
  
/* 
 * Functions can write to this file for debug info
 * Before writing to it, it should be checked whether
 * the file is not NULL:
 * if (debug) fprintf(debug,"%s","Hallo");
 */
extern FILE *debug;

extern void _range_check(int n,int n_min,int n_max,const char *var,
                         const char *file,int line);
#define range_check(n,n_min,n_max) _range_check(n,n_min,n_max,#n,__FILE__,__LINE__)
  /* Range check will terminate with an error message if not
   * n E [ n_min, n_max >
   * That is n_min is inclusive but not n_max.
   */


extern void _gmx_error(const char *key,const char *msg,const char *file,int line);
#define gmx_error(key,msg) _gmx_error(key,msg,__FILE__,__LINE__)
#define gmx_open(fn)    gmx_error("open",fn) 
#define gmx_file(msg)   gmx_error("file",msg)

#ifdef __cplusplus
	   }
#endif

#endif	/* _fatal_h */
