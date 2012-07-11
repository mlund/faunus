#ifdef ENABLE_MPI
#include <faunus/faunus.h>
#include <faunus/mpi.h>

using namespace Faunus;
using namespace std;

int main() {

  Faunus::MPI::ParticleTransmitter pt;
  Faunus::MPI::MPIController mpi;

  p_vec p(2);
  p_vec trial=p;

  if (mpi.isMaster()) {
    p[0].charge=1;
    trial=p;
    pt.sendExtra.push_back(2.0);
    pt.send(mpi, p, 1);
    pt.recv(mpi, 1, trial);
    pt.waitsend();
    pt.waitrecv();
  }
  else {
    p[0].charge=-1;
    pt.sendExtra.push_back(10.0);
    pt.send(mpi, p, 0);
    pt.recv(mpi, 0, trial);
    pt.waitsend();
    pt.waitrecv();
  }

  cout << mpi.rank << " " << p[0].charge << " " << trial[0].charge << " " << pt.recvExtra[0] << endl;
}
#endif
