#ifndef FAU_NOTIFICATION_H
#define FAU_NOTIFICATION_H

#include "faunus/common.h"

namespace Faunus {
  /*!
   * \brief Class for pushing short messages to user.
   * \author Mikael Lund
   * \date Malmo, November 2010
   * \todo Add linux support. gnotify?
   *
   * This class helps sending short messages to the user - for example
   * in the form of an email or a popup window.
   * \n\n
   * MacOS X/Growl:\n
   * User notification can be done using Growl. For this to work,
   * the 'growlnotify' command line tool must be installed and at
   * compile time the macro GROWLNOTIFYEXE must be set to the path
   * of the executable, usually /usr/local/bin/growlnotify. The build
   * system will automatically try to do this.
   */
  class notifyUser {
    private:
      string title;
    public:
      string hostname;                         //!< Hostname to send message to (if applicable)
      string email;                            //!< Email address to send message to (if applicable)
      int message(string,bool=false);          //!< Send a message
      notifyUser(string="FAUNUS Notification");
  };
} //namespace
#endif
