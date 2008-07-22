#include "imdwrap.h"

imdwrap::imdwrap(int n, int p) {
  port = p;
  npart = n;
  coords = new float[npart*3];
};

imdwrap::~imdwrap() {
  delete[] coords;
};

void imdwrap::wait_for_connection() {
  vmdsock_init();
  sock = vmdsock_create();
  clientsock = NULL;
  vmdsock_bind(sock, port);
  vmdsock_listen(sock);
  cout << "Waiting for IMD connection on port " << port << endl;
  while (!clientsock) {
    if (vmdsock_selread(sock, 0) > 0) {
      clientsock = vmdsock_accept(sock);
      if (imd_handshake(clientsock)) {
        clientsock = NULL;
      };
    }
  }
  sleep(1);
  if (vmdsock_selread(clientsock, 0) != 1 ||
      imd_recv_header(clientsock, &length) != IMD_GO) {
    clientsock = NULL;
  }
};

void imdwrap::send_particles(vector<particle> &p) {
  if (clientsock) {
    int j=0;
    imd_send_energies(clientsock, &energies);
    for (int i=0; i<npart; i++) {
      coords[j++]=p[i].x;
      coords[j++]=p[i].y;
      coords[j++]=p[i].z;
    };
    imd_send_fcoords(clientsock, npart, coords);
  }
};
