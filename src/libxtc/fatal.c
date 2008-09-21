/*
 * $Id: /local/xtcio/src/fatal.c 21 2006-08-05T20:58:43.457505Z svm  $
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

#include "fatal.h"
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include "typedefs.h"
#include "smalloc.h"

static char *fatal_tmp_file = NULL;
static FILE* stdlog = NULL;

void _where(const char *file,int line)
{
  static bool bFirst = TRUE;
  static int  nskip  = -1;
  static int  nwhere =  0;
  FILE *fp;
  char *temp; 
  
  if ( bFirst ) {
    if ((temp=getenv("WHERE")) != NULL)
      nskip = atoi(temp);
    bFirst = FALSE;
  } 

  if (nskip >= 0) {
    /* Skip the first n occasions, this allows to see where it goes wrong */
    if (nwhere >= nskip) {
      if (stdlog)
	fp = stdlog;
      else
	fp = stderr;
      fprintf(fp,"WHERE %d, file %s - line %d\n",nwhere,file,line);
    }
    nwhere++;
  }
}

static void bputc(char *msg,int *len,char ch)
{
  msg[(*len)++]=ch;
}

static void bputs(char *msg,int *len,const char *s,int fld)
{
  for (fld-=(int)strlen(s); fld>0; fld--) 
    bputc(msg,len,' ');
  while (*s) 
    bputc(msg,len,*(s++));
}

static int getfld(const char **p)
{
  int fld;

  fld=0;
  while (isdigit(**p)) fld=(fld*10)+((*((*p)++))-'0');
  return fld;
}

static int fatal_errno = 0;


static void clean_fatal_tmp_file(void);
void clean_fatal_tmp_file(void)
{
  if (fatal_tmp_file) {
    fprintf(stderr,"Cleaning up temporary file %s\n",fatal_tmp_file);
    remove(fatal_tmp_file);
    sfree(fatal_tmp_file);
    fatal_tmp_file = NULL;
  }
}

void gmx_fatal(int f_errno,const char *file,int line,const char *fmt,...)
{
  va_list ap;
  const char    *p;
  char    cval,*sval,msg[STRLEN];
  char    ibuf[64],ifmt[64];
  int     index,ival,fld,len;
  double  dval;
#ifdef _SPECIAL_VAR_ARG
  int     f_errno,line;
  char    *fmt,*file;
  
  va_start(ap,);
  f_errno = va_arg(ap,int);
  file    = va_arg(ap,char *);
  line    = va_arg(ap,int);
  fmt     = va_arg(ap,char *);
#else
  va_start(ap,fmt);
#endif

  clean_fatal_tmp_file();
  
  len=0;
  for (p=fmt; *p; p++) {
    if (*p!='%')
      bputc(msg,&len,*p);
    else {
      p++;
      fld=getfld(&p);
      switch(*p) {
      case 'x':
	ival=va_arg(ap,int);
	sprintf(ifmt,"0x%%%dx",fld);
	sprintf(ibuf,ifmt,(unsigned int)ival);
	for(index=0; (index<(int)strlen(ibuf)); index++)
	  bputc(msg,&len,ibuf[index]);
	break;
      case 'd':
	ival=va_arg(ap,int);
	sprintf(ifmt,"%%%dd",fld);
	sprintf(ibuf,ifmt,ival);
	for(index=0; (index<(int)strlen(ibuf)); index++)
	  bputc(msg,&len,ibuf[index]);
	break;
      case 'f':
	dval=va_arg(ap,double);
	sprintf(ifmt,"%%%df",fld);
	sprintf(ibuf,ifmt,dval);
	for(index=0; (index<(int)strlen(ibuf)); index++)
	  bputc(msg,&len,ibuf[index]);
	break;
      case 'c':
	cval=(char) va_arg(ap,int); /* char is promoted to int */
	bputc(msg,&len,cval);
	break;
      case 's':
	sval=va_arg(ap,char *);
	bputs(msg,&len,sval,fld);
	break;
      default:
	break;
      }
    }
  }
  va_end(ap);
  bputc(msg,&len,'\0');

  fatal_errno = f_errno;
  
  _gmx_error("fatal",msg,file,line);
}

static char   warn_buf[1024]   = "";


/* 
 * These files are global variables in the gromacs preprocessor
 * Every routine in a file that includes fatal.h can write to these
 * debug channels. Depending on the debuglevel used
 * 0 to 3 of these filed are redirected to /dev/null
 *
 */
FILE *debug=NULL;

static char *gmxuser = "Please report this to the mailing list (gmx-users@gromacs.org)";


void _gmx_error(const char *key,const char *msg,const char *file,int line)
{
  exit(1);
  /* TODO */
#if 0
  char buf[10240],tmpbuf[1024];
  int  cqnum;

  /* protect the audience from suggestive discussions */
  char *lines = "-------------------------------------------------------";
  
  cool_quote(tmpbuf,1023,&cqnum);
  sprintf(buf,"%s\nProgram %s, %s\n"
	  "Source code file: %s, line: %d\n\n"
	  "%s:\n%s\n%s\n\n%s\n",
	  lines,ShortProgram(),GromacsVersion(),file,line,
	  gmx_strerror(key),msg ? msg : warn_buf,lines,tmpbuf);
  
  gmx_error_handler(buf);
  #endif  
}

void _range_check(int n,int n_min,int n_max,const char *var,const char *file,int line)
{
  char buf[1024];
  
  if ((n < n_min) || (n >= n_max)) {
    if (strlen(warn_buf) > 0) {
      strcpy(buf,warn_buf);
      strcat(buf,"\n");
      warn_buf[0] = '\0';
    }
    else
      buf[0] = '\0';
    
    sprintf(buf+strlen(buf),"Variable %s has value %d. It should have been "
	    "within [ %d .. %d ]\n%s",var,n,n_min,n_max,gmxuser);
    
    _gmx_error("range",buf,file,line);
  }
}
