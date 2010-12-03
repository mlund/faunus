#ifndef FAU_MESSAGEPOOL_H
#define FAU_MESSAGEPOOL_H

#include <faunus/common.h>
#include <faunus/vectorid.h>

namespace Faunus {
  /*!
   * \brief Class for collecting messages from many different reporters (clients)
   *
   * This class can be used to collect short messages or "one-liners" from
   * different reporters identified by a client ID string. The reports are piled
   * together by a client id keyword. New clients are dynamically added.
   *
   * Examples\n
   * \code
   * messagepool msg;
   * msg.client["myid"]+="Minor warning about blabla!";
   * msg.client["herid"]+="Failed to open inputfile...";
   * msg.client["herid"].fatal=true;
   * cerr << msg.info();
   * if (msg.fatalerror()>-1)
   *   return 1; // bail out due to fatal error
   * \endcode
   *
   * \date Lund, 2010
   * \author Mikael Lund
   */
  class messagepool {
  public:
    class clientdata {
    public:
      bool fatal;                        //!< True if fatal error encountered
      vector<string> msg;                //!< List of reported messages
      clientdata & operator+=(string);   //!< Convenient addition of new messages
    };
    vectorid<messagepool::clientdata> client;         //!< Vector of client data
    int fatalerror();                    //!< Return index of first client w. fatal error. Else -1.
    string info();                       //!< Report with all messages grouped per client ID
  };
  
  extern messagepool errlog;  //!< Global scope
}// namespace
#endif
