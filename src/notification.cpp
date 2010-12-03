#include <faunus/notification.h>

namespace Faunus {

  /*!
   * \param titlename Title line for the message.
   */
  notifyUser::notifyUser(string titlename) {
    title=titlename;
  }

  /*!
   * \param message Text message to send to user.
   * \param sticky  If true, the message will stay visible until user reponds. Default is false.
   */
  int notifyUser::message(string message, bool sticky) {
#ifdef GROWLNOTIFYEXE
    std::ostringstream o;
    o << GROWLNOTIFYEXE << " -t \"" << title << "\" -m \"" << message << "\"";
    if (hostname.size()>0)
      o << " --host " << hostname;
    if (sticky==true)
      o << " --sticky";
    return system( o.str().c_str() );
#else
    return 0;
#endif
  }
}//namespace

