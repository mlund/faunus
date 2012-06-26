#ifdef ENABLE_MPI

#ifndef FAU_MPI_H
#define FAU_MPI_H

#include <faunus/common.h>
#include <faunus/point.h>
#include <mpi.h>

namespace Faunus {
  namespace MPI {
    class MPIController {
      public:
        MPIController(MPI_Comm=MPI_COMM_WORLD);
        ~MPIController(); //!< End of all MPI calls!
        MPI_Comm comm;    //!< Communicator (Default: MPI_COMM_WORLD)
        int nproc;        //!< Number of processors in communicator
        int rank;         //!< Rank of process
        int master;       //!< Rank number of the master
        bool isMaster();  //!< Test if current process is master
    };

    /*!
     * \brief This is a class for sending/receiving particle vectors over MPI.
     * \date Lund 2012
     * \author Mikael Lund
     *
     * Currently data is passed in single floating point precision (float).
     */
    class ParticleTransmitter {
      public:
        enum dataformat {XYZ, XYZQ};
        ParticleTransmitter();
        dataformat format;                            //!< Data format to send/receive - default is XYZQ
        void send(MPIController&, const p_vec&, int); //!< Send particle vector to another node 
        void recv(MPIController&, int, p_vec&);       //!< Receive particle vector from another node
        void waitsend();
        void waitrecv();

      private:
        vector<float> sendBuf, recvBuf;
        int tag;        //!< so far unused (zero)
        p_vec *dstPtr;  //!< pointer to receiving particle vector

        MPI_Request sendReq, recvReq;
        MPI_Status sendStat, recvStat;

        void pvec2buf(const p_vec&); //!< Copy source particle vector to send buffer
        void buf2pvec(p_vec&);       //!< Copy receive buffer to target particle vector
    };

  } //end of mpi namespace
}//end of faunus namespace

#endif
#endif
