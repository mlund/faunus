#ifdef ENABLE_MPI

#ifndef FAU_MPI_H
#define FAU_MPI_H

#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/slump.h>
#include <faunus/textio.h>
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
      private:
      public:
        MPIController(MPI_Comm=MPI_COMM_WORLD); //!< Constructor
        ~MPIController(); //!< End of all MPI calls!
        MPI_Comm comm;    //!< Communicator (Default: MPI_COMM_WORLD)
        int nproc;        //!< Number of processors in communicator
        int rank;         //!< Rank of process
        int master;       //!< Rank number of the master
        bool isMaster();  //!< Test if current process is master
        slump random;     //!< Random number generator for MPI calls
        string id;        //!< Unique name associated with current rank
        string prefix;    //!< Unique file prefix associated with current rank
        std::ofstream cout; //!< Redirect stdout to here for rank-based file output
    };

    /*!
     * \brief Class for transmitting floating point arrays over MPI
     */
    class FloatTransmitter {
      private:
        MPI_Request sendReq, recvReq;
        MPI_Status sendStat, recvStat;
        int tag;
      public:
        typedef double floatp;   //!< Transmission precision
        FloatTransmitter();
        vector<floatp> swapf(MPIController&, vector<floatp>&, int); //!< Swap data with another process
        void sendf(MPIController&, vector<floatp>&, int); //!< Send vector of floats
        void recvf(MPIController&, int, vector<floatp>&); //!< Receive vector of floats
        void waitsend(); //!< Wait for send to finish              
        void waitrecv(); //!< Wait for reception to finish
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
     *   floatp extra1 = 2.34, extra2 = -1.23
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
    class ParticleTransmitter : public FloatTransmitter {
      public:
        enum dataformat {XYZ, XYZQ};
        dataformat format;                              //!< Data format to send/receive - default is XYZQ
        vector<floatp> sendExtra;                      //!< Put extra data to send here.
        vector<floatp> recvExtra;                      //!< Received extra data will be stored here

        ParticleTransmitter();
        void send(MPIController&, const p_vec&, int); //!< Send particle vector to another node 
        void recv(MPIController&, int, p_vec&);       //!< Receive particle vector from another node
        void waitrecv();

      private:
        vector<floatp> sendBuf, recvBuf;
        p_vec *dstPtr;  //!< pointer to receiving particle vector
        void pvec2buf(const p_vec&); //!< Copy source particle vector to send buffer
        void buf2pvec(p_vec&);       //!< Copy receive buffer to target particle vector
    };

  } //end of mpi namespace
}//end of faunus namespace

#endif
#endif
