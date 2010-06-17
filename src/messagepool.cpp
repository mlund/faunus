#include <faunus/messagepool.h>
namespace Faunus {
  messagepool errlog; // Global scope

  /*!
   * \param s string containing message to add (preferably one line, only)
   */
  messagepool::clientdata & messagepool::clientdata::operator+=(string s) {
    msg.push_back(s);
    return *this;
  }
  
  int messagepool::fatalerror() {
    for (int i=0; i<client.size(); i++)
      if (client.at(i).fatal==true)
        return i;
    return -1;
  }

  string messagepool::info() {
    std::ostringstream o;
    o << endl << "# COLLECTED MESSAGES:" << endl;
    for (int i=0; i<client.size(); i++)
      if (client.at(i).msg.size()>0) {
        o << "#   " << client.id(i) << ":" << ((client.at(i).fatal==true) ? " (fatal error!)\n" : "\n");
        for (int m=0; m<client.at(i).msg.size(); m++)
          o << "#     " << client.at(i).msg[m] << endl;
      }
    return o.str();
  }
}//namespace
