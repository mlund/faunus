#include <faunus/textio.h>

namespace Faunus {
  /*!
   * Please direct all output to stdout here. By default this is exactly the same
   * as using std::cout but by using this alias it is possible to redirect all output as
   * needed in for example MPI code.
   */
  //std::ostream& textio::fcout = std::cout;

  /*!
   * As textio::fcout but for standard error.
   */
  //std::ostream& textio::fcerr = std::cerr;

  std::string textio::prefix = "";

}//namespace
