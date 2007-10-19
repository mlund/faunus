/*
 * This gives access to MacOSX intrinsic
 * functions while at the same time ensure
 * generic processor compatibility. The
 * inverse square root function may speed up
 * certain MC programs by a factor of 2.
 * 
 * Compile with -DMACOSX to activate PPC
 * intrinsic functions.
 *
 * M. Lund, 2005.
 *
*/
#ifdef MACOSX
  #include<ppc_intrinsics.h>
  #define isqrt(x) __frsqrtes(x) //estimate! carefull, may cause trouble!
#endif

// Intel ICC
#ifdef __INTEL_COMPILER
  #include <mathimf.h>
  #define isqrt(x) invsqrt(x)
#endif

// fall back to standard
//#ifndef isqrt(x)
//  #define isqrt(x) 1./sqrt(x)
//#endif

