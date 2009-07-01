// Link against libfaunus and libfortranfunc
#include "faunus/fortran.h"
int main() {
  // See: http://www.yolinux.com/TUTORIALS/LinuxTutorialMixingFortranAndC.html
  double x[10];
  int n=2;
  Faunus::gentab_();
  Faunus::vscoul_(x,x,x,&n);
}
