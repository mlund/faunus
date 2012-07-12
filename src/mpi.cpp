#ifdef ENABLE_MPI

#include <faunus/mpi.h>

namespace Faunus {
  namespace MPI {

    /*!
     * Besides initiating MPI, the current rank will be added to the global
     * file I/O prefix, textio::prefix which is useful for saving rank specific
     * data (parallel tempering, for example). The prefix format is
     * \li \c "prefix + mpi%r." where \c \%r is the rank number.
     */
    MPIController::MPIController(MPI_Comm c) : comm(c), master(0) {
      MPI_Init(NULL,NULL);
      MPI_Comm_size(comm, &nproc);
      MPI_Comm_rank(comm, &rank);    
      std::ostringstream o;
      o << rank;
      id = o.str();
      textio::prefix += "mpi" + id + ".";
      cout.open(textio::prefix+"stdout");
    }

    MPIController::~MPIController() {
      MPI_Finalize();
      cout.close();
    }

    bool MPIController::isMaster() {
      return (rank==master);
    }

    FloatTransmitter::FloatTransmitter() {
      tag=0;
    }

    void FloatTransmitter::sendf(MPIController &mpi, vector<floatp> &src, int dst) {
      MPI_Issend(&src[0], src.size(), MPI_DOUBLE, dst, tag, mpi.comm, &sendReq);
    }

    void FloatTransmitter::waitsend() {
      MPI_Wait(&sendReq, &sendStat);
    } 

    void FloatTransmitter::recvf(MPIController &mpi, int src, vector<floatp> &dst) {
      MPI_Irecv(&dst[0], dst.size(), MPI_DOUBLE, src, tag, mpi.comm, &recvReq);
    }

    void FloatTransmitter::waitrecv() {
      MPI_Wait(&recvReq, &recvStat);
    }

    /*!
     * This will send a vector of floats and at the same time wait for the destination process
     * to send back another vector of the same size.
     * 
     * \param mpi MPI controller to use
     * \param src Vector to send
     * \param dst Node to send/receive to/from
     *
     * \todo Use MPI_Sendrecv( sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, &status);
     */
    vector<FloatTransmitter::floatp> FloatTransmitter::swapf(MPIController &mpi, vector<floatp> &src, int dst) {
      vector<floatp> v( src.size() );
      recvf(mpi, dst, v);
      sendf(mpi, src, dst);
      waitrecv();
      waitsend();
      return v;
    }

    ParticleTransmitter::ParticleTransmitter() {
      format=XYZQ;
    }


    /*!
     * \param mpi MPI controller to use
     * \param src Source particle vector
     * \param dst Destination node
     */
    void ParticleTransmitter::send(MPIController &mpi, const p_vec &src, int dst) {
      assert(dst>=0 && dst<mpi.nproc && "Invalid MPI destination");
      pvec2buf(src);
      FloatTransmitter::sendf(mpi, sendBuf, dst);
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

    /*!
     * \param mpi MPI controller to use
     * \param src Source node
     * \param dst Destination particle vector
     */
    void ParticleTransmitter::recv(MPIController &mpi, int src, p_vec &dst) {
      assert(src>=0 && src<mpi.nproc && "Invalid MPI source");
      dstPtr=&dst;   // save a pointer to the destination particle vector
      if (format==XYZ)
        recvBuf.resize( 3*dst.size() );
      if (format==XYZQ)
        recvBuf.resize( 4*dst.size() );

      // resize to fit extra data (if any)
      recvExtra.resize( sendExtra.size() );
      recvBuf.resize( recvBuf.size() + recvExtra.size() );

      FloatTransmitter::recvf(mpi, src, recvBuf);
    }

    void ParticleTransmitter::waitrecv() {
      FloatTransmitter::waitrecv();
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
      if ( (size_t)i!=recvBuf.size() )
        cout << "!!!!!!!!!!!\n";
    }

  } //end of mpi namespace
}//end of faunus namespace

#endif
