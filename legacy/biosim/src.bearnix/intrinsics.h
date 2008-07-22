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
#ifdef XLC
//IBM XL C++ yet to be implemented...
#endif

#ifdef MACOSX
 #include<ppc_intrinsics.h>
 #define frsqrte(x) __frsqrtes(x) //estimate! carefull, may cause trouble!
// #define pow2(x) __fmuls(x,x)
// #define fmul(x,y) __fmuls(x,y)
#else
 #define frsqrte(x) 1./sqrt(x)
// #define pow2(x) x*x
// #define fmul(x,y) x*y
#endif

