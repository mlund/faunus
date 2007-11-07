/*
 */

#ifndef _imdwrap_h
#define _imdwrap_h

#include "vmd/vmdsock.h"
#include "vmd/imd.h"
#include "point.h"
#include <vector>
#include <iostream>

using namespace std;

class imdwrap {
  private:
    int port;  // imd port number
    int  npart; // number of particles
    void *sock;
    void *clientsock;
    int length;
    IMDEnergies energies;
    float * coords;
  public:
    imdwrap(int, int=54321);
    ~imdwrap();
    void wait_for_connection();
    void send_particles(vector<particle> &);
};

#endif
