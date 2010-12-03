#ifndef FAU_COUNTDOWN_H
#define FAU_COUNTDOWN_H

#ifndef SWIG
#include <string>
#include <ctime>
#endif

namespace Faunus {

/*! \brief Estimate speed of a computational process
 *  \author Mikael Lund
 * 
 *  This class can be used to estimate when a computational
 *  process will finish. The current resolution is seconds.
 */
template<class T>
class countdown {
  private:
    T max;
    int time_i;              //!< Starting time in seconds
    time_t rawtime; 
    struct tm *timeinfo;
  public:
    countdown(T);
    float speed(T);          //!< Calculate speed
    int elapsed();           //!< Time elapsed in seconds
    std::string eta(T);      //!< Estimate time of arrival
};
/*! \param maxvalue Value at arrival
 */
template<class T> countdown<T>::countdown(T maxvalue) {
  time_i=time(0);
  max=maxvalue;
}
template<class T> float countdown<T>::speed(T midvalue) {
  return float(time(0)-time_i) / midvalue;
}
/*! \param midvalue Value somewhere between start and arrival
 *  \return String with estimated time and date of arrival
 */
template<class T> std::string countdown<T>::eta(T midvalue) {
  rawtime = time(NULL) + int( speed(midvalue) * (max-midvalue)  );
  timeinfo = localtime(&rawtime);
  return asctime(timeinfo);
}
template<class T> int countdown<T>::elapsed() { return time(0)-time_i; }

}//namespace
#endif
