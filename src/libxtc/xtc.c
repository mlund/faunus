/*
 * $Id: /local/xtcio/src/xtc.c 39 2008-09-14T06:31:00.970588Z lidb  $
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

#include "xtc.h"

#include <assert.h>

#include "smalloc.h"
#include "gmxfio.h"
#include "fatal.h"

#define DIM 3
#define XX 0
#define YY 1
#define ZZ 2

struct Xtc {
  XDR m_xdr;
  FILE* m_file;
};
static const int header_size = 16;

#ifdef HAVE_FSEEKO
#define gmx_fseek(A,B,C) fseeko(A,B,C)
#define gmx_ftell(A) ftello(A)
#define gmx_off_t off_t
#else
#define gmx_fseek(A,B,C) fseek(A,B,C)
#define gmx_ftell(A) ftell(A)
#define gmx_off_t int
#endif

static int xtc_check(char *str,bool bResult,char *file,int line);
#define XTC_CHECK(s,b) xtc_check(s,b,__FILE__,__LINE__)

static int xdr_real(XDR *xdrs,
                    xtc_real *r)
{
#ifdef GMX_DOUBLE
    float f;
    int   ret;
    
    f=*r;
    ret=xdr_float(xdrs,&f);
    *r=f;
    
    return ret;
#else
    return xdr_float(xdrs,(float *)r);
#endif
}

static void xtc_init(Xtc* self) {
  self->m_file = 0;
}

static const char*
get_newmode(const char* mode) {
  switch(*mode) {
    case 'w':
    case 'W':
      return "wb+";
    case 'a':
    case 'A':
      return "ab+";
    case 'r':
    case 'R':
      return "rb";
    default:
      return NULL;
  }
}

static enum xdr_op
get_op(const char* mode) {
  switch(*mode) {
    case 'w':
    case 'W':
    case 'a':
    case 'A':
      return XDR_ENCODE;
    case 'r':
    case 'R':
      return XDR_DECODE;
    default:
      return (enum xdr_op) -1;
  }
}

Xtc* xtc_create(const char* filename, const char* mode) {
  Xtc* res = 0;
  const char* newtype = get_newmode(mode);
  enum xdr_op op = get_op(mode);
  FILE* file = 0;
  
  if(newtype == NULL) 
    return NULL;
  assert(op == XDR_ENCODE
         || op == XDR_DECODE);
    
  file = fopen(filename, newtype); 
  if(file == NULL) {
    return NULL;
  }

  snew(res, 1);
  xtc_init(res);
  res->m_file = file;
  xdrstdio_create(&res->m_xdr,
                  res->m_file,
                  op);

  return res;
}

Xtc* xtc_stdio_create(FILE* file, const char* mode)
{
  enum xdr_op op = get_op(mode);
  if(op != XDR_ENCODE
      && op != XDR_DECODE) return NULL;
  
  Xtc* res = 0;
  xtc_init(res);
  res->m_file = file;
  xdrstdio_create(&res->m_xdr,
                  res->m_file,
                  op);
  return res;
}

void xtc_destroy(Xtc* self)
{
  xdr_destroy(&self->m_xdr);
  if(self->m_file) {
    fclose(self->m_file);
  }
}


static void check_xtc_magic(int magic)
{
  if (magic != XTC_MAGIC) 
    gmx_fatal(FARGS,"Magic Number Error in XTC file (read %d, should be %d)",
		magic,XTC_MAGIC);
}

int xtc_check(char *str,bool bResult,char *file,int line)
{
  if (!bResult) {
    if (debug)
      fprintf(debug,"\nXTC error: read/write of %s failed, "
	      "source file %s, line %d\n",str,file,line);
    return 0;
  }
  return 1;
}


static int xtc_header(XDR *xd,int *magic,int *natoms,int *step,xtc_real *time,
		      bool *bOK)
{
  int result;

  if (xdr_int(xd,magic) == 0)
    return 0;
  result=XTC_CHECK("natoms", xdr_int(xd,natoms));  /* number of atoms */
  if (result)
    result=XTC_CHECK("step",   xdr_int(xd,step));    /* frame number    */
  if (result)
    result=XTC_CHECK("time",   xdr_real(xd,time));   /* time            */
  *bOK=(result!=0);

  return result;
}

static int xtc_coord(XDR *xd,
                     int *natoms,
                     xtc_matrix box,
                     xtc_rvec *x,
                     xtc_real *prec,
                     bool bRead)
{
  int i,j,result;
#ifdef GMX_DOUBLE
  float *ftmp;
  float fprec;
#endif
    
  /* box */
  result=1;
  for(i=0; ((i<DIM) && result); i++)
    for(j=0; ((j<DIM) && result); j++)
      result=XTC_CHECK("box",xdr_real(xd,&(box[i][j])));
  
  if (!result)
      return result;
  
#ifdef GMX_DOUBLE
  /* allocate temp. single-precision array */
  snew(ftmp,(*natoms)*DIM);
  
  /* Copy data to temp. array if writing */
  if(!bRead)
  {
      for(i=0; (i<*natoms); i++)
      {
          ftmp[DIM*i+XX]=x[i][XX];      
          ftmp[DIM*i+YY]=x[i][YY];      
          ftmp[DIM*i+ZZ]=x[i][ZZ];      
      }
      fprec = *prec;
  }
  result=XTC_CHECK("x",xdr3dfcoord(xd,ftmp,natoms,&fprec));
  
  /* Copy from temp. array if reading */
  if(bRead)
  {
      for(i=0; (i<*natoms); i++)
      {
          x[i][XX] = ftmp[DIM*i+XX];      
          x[i][YY] = ftmp[DIM*i+YY];      
          x[i][ZZ] = ftmp[DIM*i+ZZ];      
      }
      *prec = fprec;
  }  
  sfree(ftmp);
#else
    result=XTC_CHECK("x",xdr3dfcoord(xd,x[0],natoms,prec)); 
#endif 
    
  return result;
}



bool xtc_write(Xtc* self,
              int natoms,
              int step,
              xtc_real time,
              xtc_matrix box,
              xtc_rvec *x,
              xtc_real prec)
{
  int magic_number = XTC_MAGIC;
  XDR *xd;
  bool bDum;
  int fp = fileno(self->m_file);

  xd = fio_getxdr(fp);
  /* write magic number and xtc identidier */
  if (!xtc_header(xd,&magic_number,&natoms,&step,&time,&bDum))
    return 0;
    
  /* write data */
  return xtc_coord(xd,&natoms,box,x,&prec,FALSE);
}

bool xtc_read_first(Xtc* self,
                   int *natoms,
                   int *step,
                   xtc_real *time,
                   xtc_matrix box,
                   xtc_rvec **x,
                   xtc_real *prec)
/*
int xtc_read_first(int fp,int *natoms,int *step,real *time,
matrix box,rvec **x,real *prec,bool *bOK) */
{
  int magic;
  XDR *xd;
  
  bool bOK=TRUE;
  int fp = fileno(self->m_file);
  xd = fio_getxdr(fp);
  
  /* read header and malloc x */
  if ( !xtc_header(xd,&magic,natoms,step,time,&bOK))
    return 0;
    
  /* Check magic number */
  check_xtc_magic(magic);
  
  snew(*x,*natoms);

  bOK=xtc_coord(xd,natoms,box,*x,prec,TRUE);
  
  return bOK;
}

bool xtc_read_next(Xtc* self,
                  int natoms,
                  int *step,
                  xtc_real *time,
                  xtc_matrix box,
                  xtc_rvec *x,
                  xtc_real *prec)
{
  int magic;
  int n;
  XDR *xd;

  bool bOK=TRUE;
  int fp = fileno(self->m_file);

  xd = fio_getxdr(fp);
  
  /* read header */
  if (!xtc_header(xd,&magic,&n,step,time,&bOK))
    return 0;
  if (n>natoms)
    gmx_fatal(FARGS, "Frame contains more atoms (%d) than expected (%d)", 
		n, natoms);
    
  /* Check magic number */
  check_xtc_magic(magic);

  bOK=xtc_coord(xd,&natoms,box,x,prec,TRUE);

  return bOK;
}


double
xtc_get_last_frame_time(Xtc* self,
                        int natoms, 
                        bool* bOK)
{
    float  time;
    gmx_off_t  off;
    int res;
    int fp = fileno(self->m_file);

    *bOK = 1;
    off = gmx_ftell(xdrfiles[fp+1]);  
    if(off < 0){
      *bOK = 0;
      return -1;
    }
    
    if( (res = gmx_fseek(xdrfiles[fp+1],-4,SEEK_END)) != 0){
      *bOK = 0;
      return -1;
    }

    time = _xtc_get_current_frame_time(fp, natoms, bOK);
    if(!(*bOK)){
      return -1;
    }
    
    if( (res = gmx_fseek(xdrfiles[fp+1],off,SEEK_SET)) != 0){
      *bOK = 0;
      return -1;
    } 
    return time;
}


int xtc_get_last_frame_number(Xtc* self, 
                          int natoms, 
                          bool* bOK)
{
    int    frame;
    gmx_off_t  off;
    int fp = fileno(self->m_file);

    *bOK = 1;
    
    if((off = gmx_ftell(xdrfiles[fp+1])) < 0){
      *bOK = 0;
      return -1;
    }

    
    if(gmx_fseek(xdrfiles[fp+1],-4,SEEK_END)){
      *bOK = 0;
      return -1;
    }

    frame = _xtc_get_current_frame_number(fp, natoms, bOK);
    if(!bOK){
      return -1;
    }


    if(gmx_fseek(xdrfiles[fp+1],off,SEEK_SET)){
      *bOK = 0;
      return -1;
    }    

    return frame;
}

bool
xtc_seek_frame(Xtc* self,
               int frame, 
               int natoms)
{
    gmx_off_t low = 0;
    gmx_off_t high,pos;
    int fp = fileno(self->m_file);

    
    /* round to 4 bytes */
    int fr;
    gmx_off_t  offset;
    if(gmx_fseek(xdrfiles[fp+1],0,SEEK_END)){
      return false;
    }

    if((high = gmx_ftell(xdrfiles[fp+1])) < 0){
      return false;
    }
    
    /* round to 4 bytes  */
    high /= sizeof(int);
    high *= sizeof(int);
    offset = ((high/2)/sizeof(int))*sizeof(int);
    
    if(gmx_fseek(xdrfiles[fp+1],offset,SEEK_SET)){
      return false;
    }
    
    while(1)
    {
        fr = _xtc_get_next_frame_number(fp,natoms);
        if(fr < 0)
        {
            return false;
        }
        if(fr != frame && abs(low-high) > header_size)
        {
            if(fr < frame)
            {
                low = offset;      
            }
            else
            {
                high = offset;      
            }
            /* round to 4 bytes */
            offset = (((high+low)/2)/4)*4;
            
            if(gmx_fseek(xdrfiles[fp+1],offset,SEEK_SET)){
	      return false;
	    }
        }
        else
        {
            break;
        }
    }
    if(offset <= header_size)
    {
        offset = low;
    }
    
    if(gmx_fseek(xdrfiles[fp+1],offset,SEEK_SET)){
      return false;
    }

    if((pos = _xtc_get_next_frame_start(fp,natoms))< 0){
    /* we probably hit an end of file */
      return false;
    }
    
    if(gmx_fseek(xdrfiles[fp+1],pos,SEEK_SET)){
      return false;
    }
    
    return true;
}

     

bool xtc_seek_time(Xtc* self, 
              double time, 
              int natoms)
{
    float t;
    float dt;
    bool bOK;
    gmx_off_t low = 0;
    gmx_off_t high,offset,pos;
    int dt_sign = 0;
    int fp = fileno(self->m_file);

  
    if(gmx_fseek(xdrfiles[fp+1],0,SEEK_END)){
      return false;
    }
    
    if((high = gmx_ftell(xdrfiles[fp+1])) < 0){
      return false;
    }
    /* round to int  */
    high /= sizeof(int);
    high *= sizeof(int);
    offset = ((high/2)/sizeof(int))*sizeof(int);

    if(gmx_fseek(xdrfiles[fp+1],offset,SEEK_SET)){
      return false;	
   }

    dt = _xtc_estimate_dt(fp,natoms,&bOK);
      
      
    
    if(!bOK)
    {
        return false;
    }
    
    while(1)
    {
	dt = _xtc_estimate_dt(fp,natoms,&bOK);
        if(!bOK)
        {
            return false;
        }else{
	  if(dt > 0){
	    if(dt_sign == -1){
	      /* Found a place in the trajectory that has positive time step while
		 other has negative time step */
	      return false;
	    }
	    dt_sign = 1;
	  }else if(dt < 0){
	    if(dt_sign == 1){
	      /* Found a place in the trajectory that has positive time step while
		 other has negative time step */
	      return false;
	    }
	    dt_sign = -1;
	  }	  
	}
        t = _xtc_get_next_frame_time(fp,natoms,&bOK);
        if(!bOK)
        {
            return false;
        }

	/* If we are before the target time and the time step is positive or 0, or we have
	 after the target time and the time step is negative, or the difference between 
	the current time and the target time is bigger than dt and above all the distance between high
	and low is bigger than 1 frame, then do another step of binary search. Otherwise stop and check
	if we reached the solution */
        if((((t < time && dt_sign >= 0) || (t > time && dt_sign == -1)) || ((t-time) >= dt && dt_sign >= 0) || ((time-t) >= -dt && dt_sign < 0))  && (abs(low-high) > header_size))
        {
	  if(dt >= 0 && dt_sign != -1)
            {
                if(t < time)
                {
                    low = offset;      
                }
                else
                {
                    high = offset;      
                }
	    }
	  else if(dt <= 0 && dt_sign == -1)
            {
	      if(t >= time)
                {
		  low = offset;      
                }
	      else
                {
		  high = offset;      
                }
	    }else{
	      /* We should never reach here */
	      return false;
	    }
            /* round to 4 bytes and subtract header*/
            offset = (((high+low)/2)/sizeof(int))*sizeof(int);
            if(gmx_fseek(xdrfiles[fp+1],offset,SEEK_SET)){
	      return false;
	    }
        }
        else
        {
            if(abs(low-high) <= header_size)
            {
                break;
            }
            /* reestimate dt */
            if(_xtc_estimate_dt(fp,natoms,&bOK) != dt)
            {
                if(bOK)
                {
                    dt = _xtc_estimate_dt(fp,natoms,&bOK);
                }
            }
            if(t >= time && t-time < dt)
            {
                break;
            }
        }        
    }
    
    if(offset <= header_size)
    {
        offset = low;
    }
    
    gmx_fseek(xdrfiles[fp+1],offset,SEEK_SET);

    if((pos = _xtc_get_next_frame_start(fp,natoms)) < 0){
      return false;
    }
    
    if(gmx_fseek(xdrfiles[fp+1],pos,SEEK_SET)){
      return false;
    }
    return true;
}
