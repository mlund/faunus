#ifdef ENABLE_MPI

#ifndef FAU_MPI_H
#define FAU_MPI_H

#include <faunus/common.h>
#include <faunus/point.h>
#include <mpi.h>

namespace Faunus {

  /*!
   * \brief Namespace for Message Parsing Interface (MPI) functionality
   */
  namespace MPI {

    /*!
     * \brief Main controller for MPI calls
     */
    class MPIController {
      public:
        MPIController(MPI_Comm=MPI_COMM_WORLD); //!< Constructor
        ~MPIController(); //!< End of all MPI calls!
        MPI_Comm comm;    //!< Communicator (Default: MPI_COMM_WORLD)
        int nproc;        //!< Number of processors in communicator
        int rank;         //!< Rank of process
        int master;       //!< Rank number of the master
        bool isMaster();  //!< Test if current process is master
    };

    /*!
     * \brief Class for sending/receiving particle vectors over MPI.
     * \date Lund 2012
     * \author Mikael Lund
     * \note Data is passed in single floating point precision (MPI_FLOAT).
     *
     * This will take a particle vector and send selected information though MPI. It is
     * possible to send only coordinates using the dataformat XYZ or, if charges should be
     * send too, XYZQ.
     *
     * Besides particle data it is possible to send extra floats by adding
     * these to the sendExtra vector; received extras will be stored in recvExtra. Before
     * transmitting extra data, make sure that recvExtra and sendExtra have the
     * same size.
     *
     * \code
     *   int dst_rank = 1;
     *   float extra1 = 2.34, extra2 = -1.23
     *   p_vec myparticles(200); // we have 200 particles
     *
     *   Faunus::MPI::MPIController mpi;
     *   Faunus::MPI::ParticleTransmitter pt;
     *
     *   pt.sendExtra.push_back(extra1);
     *   pt.sendExtra.push_back(extra2);
     *
     *   pt.send(mpi, myparticles, dst_rank);
     *   pt.waitsend();
     * \endcode
     */
    class ParticleTransmitter {
      public:
        enum dataformat {XYZ, XYZQ};
        dataformat format;                            //!< Data format to send/receive - default is XYZQ
        vector<float> sendExtra;                      //!< Put extra data to send here.
        vector<float> recvExtra;                      //!< Received extra data will be stored here

        ParticleTransmitter();
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
