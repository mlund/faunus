#ifndef _stopwatch_h
#define _stopwatch_h

#include <string>

/*! \brief Estimate speed of a computational process
 *  \author Mikael Lund
 * 
 *  This class can be used to estimate when a computational
 *  process will finish. The current resolution is seconds.
 */
template<class T>
class stopwatch {
  private:
    T max;
    int time_i;
    time_t rawtime; 
    struct tm *timeinfo;
  public:
    stopwatch(T);
    float speed(T);          //!< Calculate speed
    std::string eta(T);      //!< Estimate time of arrival
};

/*! \param maxvalue Value at arrival
 */
template<class T>
stopwatch<T>::stopwatch(T maxvalue) {
  time_i=time(0);
  max=maxvalue;
}

template<class T>
float stopwatch<T>::speed(T midvalue) {
  return float(time(0)-time_i) / midvalue;
}

/*! \param midvalue Value somewhere between start and arrival
 *  \return String with estimated time and date of arrival
 */
template<class T>
std::string stopwatch<T>::eta(T midvalue) {
  rawtime = time(NULL) + int( speed(midvalue) * (max-midvalue)  );
  timeinfo = localtime(&rawtime);
  return asctime(timeinfo);
}

#endif
