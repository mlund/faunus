#ifdef ENABLE_MPI

#include <faunus/mpi.h>

namespace Faunus {
  namespace MPI {

    MPIController::MPIController(MPI_Comm c) : comm(c), master(0) {
      MPI_Init(NULL,NULL);
      MPI_Comm_size(comm, &nproc);
      MPI_Comm_rank(comm, &rank);    
    }

    MPIController::~MPIController() {
      MPI_Finalize();
    }

    bool MPIController::isMaster() {
      return (rank==master);
    }

    ParticleTransmitter::ParticleTransmitter() {
      tag=0;
      format=XYZQ;
    }

    /*!
     * \param mpi MPI controller to use
     * \param src Source particle vector
     * \param dst Destination node
     */
    void ParticleTransmitter::send(MPIController &mpi, const p_vec &src, int dst) {
      assert(dst<mpi.nproc && "Invalid MPI destination");
      pvec2buf(src);
      MPI_Issend(&sendBuf[0], sendBuf.size(), MPI_FLOAT, dst, tag, mpi.comm, &sendReq);
    }

    void ParticleTransmitter::pvec2buf(const p_vec &src) {
      sendBuf.clear();
      for (auto &p : src) {
        sendBuf.push_back(p.x);
        sendBuf.push_back(p.y);
        sendBuf.push_back(p.z);
        if (format==XYZQ)
          sendBuf.push_back(p.charge);
      }
      for (auto i : sendExtra)
        sendBuf.push_back(i);
    }

    void ParticleTransmitter::waitsend() {
      MPI_Wait(&sendReq, &sendStat);
    }

    /*!
     * \param mpi MPI controller to use
     * \param src Source node
     * \param dst Destination particle vector
     */
    void ParticleTransmitter::recv(MPIController &mpi, int src, p_vec &dst) {
      assert(src<mpi.nproc && "Invalid MPI source");
      dstPtr=&dst;   // save a pointer to the destination particle vector
      if (format==XYZ)
        recvBuf.resize( 3*dst.size() );
      if (format==XYZQ)
        recvBuf.resize( 4*dst.size() );

      // resize to fit extra data (if any)
      recvExtra.resize( sendExtra.size() );
      recvBuf.resize( recvBuf.size() + recvExtra.size() );

      MPI_Irecv(&recvBuf[0], recvBuf.size(), MPI_FLOAT, src, tag, mpi.comm, &recvReq);
    }

    void ParticleTransmitter::waitrecv() {
      MPI_Wait(&recvReq, &recvStat);
      buf2pvec(*dstPtr);
    }

    void ParticleTransmitter::buf2pvec(p_vec &dst) {
      int i=0;
      for (auto &p : dst) {
        p.x=recvBuf[i++];
        p.y=recvBuf[i++];
        p.z=recvBuf[i++];
        if (format==XYZQ)
          p.charge=recvBuf[i++];
      }
      for (auto &x : recvExtra)
        x=recvBuf[i++];
      assert( (size_t)i==recvBuf.size() );
    }

  } //end of mpi namespace
}//end of faunus namespace

#endif
