#ifndef faunus_fortran_h
#define faunus_fortran_h
// Definitions for the fortran routines in libfortranfunc library
// (see legacy/fortran)
namespace Faunus {
  extern "C"
  {
    void vscoul_(double*, double*, double*, int* );
    void gentab_();
  }
};
#endif
