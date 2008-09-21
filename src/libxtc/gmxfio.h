/*
 * $Id: /local/xtcio/src/gmxfio.h 21 2006-08-05T20:58:43.457505Z svm  $
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

#ifndef _gmxfio_h
#define _gmxfio_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "xdrf.h"

/********************************************************
 * Open and Close 
 ********************************************************/

extern int fio_open(const char *fn, const char *mode);
/* Open a new file for reading or writing.
 * The file type will be deduced from the file name.
 * If fn is NULL, stdin / stdout will be used for Ascii I/O (TPA type)
 * mode may be "r", "w", or "a". You should append a "b" to the mode
 * if you are writing a binary file, but the routine will also 
 * doublecheck it and try to do it if you forgot. This has no effect on
 * unix, but is important on windows.
 */
 
extern void fio_close(int fp);
/* Close the file corresponding to fp (if not stdio)
 * The routine will exit when an invalid fio is handled.
 */


extern XDR *fio_getxdr(int fio);
/* Return the file pointer itself */

#endif
