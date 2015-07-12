#ifdef ENABLE_MPI

#ifndef FAU_MPI_H
#define FAU_MPI_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/slump.h>
#include <faunus/textio.h>
#include <mpi.h>
#endif

namespace Faunus {

  /**
   * @brief Namespace for Message Parsing Interface (MPI) functionality
   */
  namespace MPI {

    /**
     * @brief Main controller for MPI calls
     *
     * This is the MPI controller required for all MPI programs.
     *
     *     MPIController mpi; // call this very first thing in your program
     *     std::cout << "I'm rank " << mpi.rank << " out of " << mpi.nproc;
     *     mpi.cout << "This will go to a file called mpi%r.stdout where %r is my rank"
     *     if (mpi.isMaster())
     *       cout << "I'm the master!";
     *
     * When MPIController is instantiated the textio::prefix variable is automatically
     * set to `mpi%j`. Which can be used to prefix input and output files. For example:
     *
     *     InputMap mcp(textio::prefix+"input"); // tries to load "mpi%r.input" where %r is the rank
     *
     * @date Lund 2012
     */
    class MPIController {
      public:
        MPIController(MPI_Comm=MPI_COMM_WORLD); //!< Constructor
        ~MPIController(); //!< End of all MPI calls!
        MPI_Comm comm;    //!< Communicator (Default: MPI_COMM_WORLD)
        int nproc();      //!< Number of processors in communicator
        int rank();       //!< Rank of process
        int rankMaster(); //!< Rank number of the master
        bool isMaster();  //!< Test if current process is master
        RandomTwister<> random; //!< Random number generator for MPI calls
        string id;        //!< Unique name associated with current rank
        std::ofstream cout; //!< Redirect stdout to here for rank-based file output

        inline string info() {
          std::ostringstream o;
          o << textio::header("Message Parsing Interface (MPI)")
            << textio::pad(textio::SUB, 25, "Number of processors") << nproc() << endl
            << textio::pad(textio::SUB, 25, "Current rank") << rank() << endl
            << textio::pad(textio::SUB, 25, "Master rank") << rankMaster() << endl;
          return o.str();
        }

      private:
        int _nproc;        //!< Number of processors in communicator
        int _rank;         //!< Rank of process
        int _master;       //!< Rank number of the master
    };

    /**
     * @brief Split N items into nproc parts
     *
     * This returns a pair with the first and last
     * item for the current rank.
     */
    template<class T=int>
      std::pair<T,T> splitEven(MPIController &mpi, T N) {
        T M = mpi.nproc();
        T i = mpi.rank();
        T beg=(N*i)/M;
        T end=(N*i+N)/M-1;
        return std::pair<T,T>(beg,end);
      }

    /**
     * @brief Reduced sum
     *
     * Each rank sends "local" to master who sums them up.
     * Master sends back (broadcasts) sum to all ranks.
     */
    inline double reduceDouble(MPIController &mpi, double local) {
      double sum;
      MPI_Allreduce(&local,&sum,1,MPI_DOUBLE,MPI_SUM,mpi.comm);
      return sum;
    }

    /*!
     * \brief Class for transmitting floating point arrays over MPI
     * \note If you change the floatp typedef, remember also to change to change to/from
     *       MPI_FLOAT or MPI_DOUBLE.
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

    /**
     * @brief Class for sending/receiving particle vectors over MPI.
     *
     * This will take a particle vector and send selected information though MPI. It is
     * possible to send only coordinates using the dataformat `XYZ` or, if charges should be
     * send too, `XYZQ`.
     *
     * Besides particle data it is possible to send extra floats by adding
     * these to the `sendExtra` vector; received extras will be stored in `recvExtra`. Before
     * transmitting extra data, make sure that `recvExtra` and `sendExtra` have the
     * same size.
     *
     *     int dst_rank = 1;
     *     floatp extra1 = 2.34, extra2 = -1.23
     *     p_vec myparticles(200); // we have 200 particles
     *    
     *     Faunus::MPI::MPIController mpi;
     *     Faunus::MPI::ParticleTransmitter pt;
     *    
     *     pt.sendExtra.push_back(extra1);
     *     pt.sendExtra.push_back(extra2);
     *    
     *     pt.send(mpi, myparticles, dst_rank);
     *     pt.waitsend();
     *
     * @date Lund 2012
     *
     */
    template<typename p_vec>
      class ParticleTransmitter : public FloatTransmitter {
        private:
          typedef typename p_vec::value_type particle; 
        public:
          enum dataformat {XYZ=3, XYZQ=4, XYZQI=5};
          vector<floatp> sendExtra;                      //!< Put extra data to send here.
          vector<floatp> recvExtra;                      //!< Received extra data will be stored here

          ParticleTransmitter();
          void send(MPIController&, const p_vec&, int); //!< Send particle vector to another node 
          void recv(MPIController&, int, p_vec&);       //!< Receive particle vector from another node
          void waitrecv();
          void setFormat(dataformat);
          void setFormat(string);
          dataformat getFormat();

        private:
          dataformat format;                             //!< Data format to send/receive - default is XYZQ
          vector<floatp> sendBuf, recvBuf;
          p_vec *dstPtr;  //!< pointer to receiving particle vector
          void pvec2buf(const p_vec&); //!< Copy source particle vector to send buffer
          void buf2pvec(p_vec&);       //!< Copy receive buffer to target particle vector
      };

    /*!
     * Besides initiating MPI, the current rank will be added to the global
     * file I/O prefix, textio::prefix which is useful for saving rank specific
     * data (parallel tempering, for example). The prefix format is
     * \li \c "prefix + mpi%r." where \c \%r is the rank number.
     */
    MPIController::MPIController(MPI_Comm c) : comm(c), _master(0) {
      MPI_Init(NULL,NULL);
      MPI_Comm_size(comm, &_nproc);
      MPI_Comm_rank(comm, &_rank);
      id=std::to_string(_rank);
      textio::prefix += "mpi" + id + ".";
      cout.open(textio::prefix+"stdout");
    }

    MPIController::~MPIController() {
      MPI_Finalize();
      cout.close();
    }

    int MPIController::nproc() { return _nproc; }

    int MPIController::rank() { return _rank; }

    int MPIController::rankMaster() { return _master; }

    bool MPIController::isMaster() { return (_rank==_master); }

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

    /**
     * This will send a vector of floats and at the same time wait for the destination process
     * to send back another vector of the same size.
     * 
     * @param mpi MPI controller to use
     * @param src Vector to send
     * @param dst Node to send/receive to/from
     *
     * @todo Use MPI_Sendrecv( sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount,
     *       recvtype, source, recvtag, comm, &status);
     */
    vector<FloatTransmitter::floatp> FloatTransmitter::swapf(MPIController &mpi, vector<floatp> &src, int dst) {
      vector<floatp> v( src.size() );
      recvf(mpi, dst, v);
      sendf(mpi, src, dst);
      waitrecv();
      waitsend();
      return v;
    }

    template<typename p_vec>
      ParticleTransmitter<p_vec>::ParticleTransmitter() { setFormat(XYZQI); }

    template<typename p_vec>
      void ParticleTransmitter<p_vec>::setFormat(dataformat d) { format = d; }

    template<typename p_vec>
      void ParticleTransmitter<p_vec>::setFormat(string s) {
        setFormat(XYZQI);
        if (s=="XYZQ")
          setFormat(XYZQ);
        if (s=="XYZ")
          setFormat(XYZ);
      }

    template<typename p_vec>
      typename ParticleTransmitter<p_vec>::dataformat
      ParticleTransmitter<p_vec>::getFormat() { return format; }

    /*!
     * \param mpi MPI controller to use
     * \param src Source particle vector
     * \param dst Destination node
     */
    template<typename p_vec>
      void ParticleTransmitter<p_vec>::send(MPIController &mpi, const p_vec &src, int dst) {
        assert(dst>=0 && dst<mpi.nproc() && "Invalid MPI destination");
        pvec2buf(src);
        FloatTransmitter::sendf(mpi, sendBuf, dst);
      }

    template<typename p_vec>
      void ParticleTransmitter<p_vec>::pvec2buf(const p_vec &src) {
        sendBuf.clear();
        for (auto &p : src) {
          sendBuf.push_back(p.x());
          sendBuf.push_back(p.y());
          sendBuf.push_back(p.z());
          if (format==XYZQ)
            sendBuf.push_back(p.charge);
          if (format==XYZQI) {
            sendBuf.push_back(p.charge);
            sendBuf.push_back( (floatp)p.id );
          }
        }
        for (auto i : sendExtra)
          sendBuf.push_back(i);
      }

    /*!
     * \param mpi MPI controller to use
     * \param src Source node
     * \param dst Destination particle vector
     */
    template<typename p_vec>
      void ParticleTransmitter<p_vec>::recv(MPIController &mpi, int src, p_vec &dst) {
        assert(src>=0 && src<mpi.nproc() && "Invalid MPI source");
        dstPtr=&dst;   // save a pointer to the destination particle vector
        if (format==XYZ)
          recvBuf.resize( 3*dst.size() );
        if (format==XYZQ)
          recvBuf.resize( 4*dst.size() );
        if (format==XYZQI)
          recvBuf.resize( 5*dst.size() );

        // resize to fit extra data (if any)
        recvExtra.resize( sendExtra.size() );
        recvBuf.resize( recvBuf.size() + recvExtra.size() );

        FloatTransmitter::recvf(mpi, src, recvBuf);
      }

    template<typename p_vec>
      void ParticleTransmitter<p_vec>::waitrecv() {
        FloatTransmitter::waitrecv();
        buf2pvec(*dstPtr);
      }

    template<typename p_vec>
      void ParticleTransmitter<p_vec>::buf2pvec(p_vec &dst) {
        int i=0;
        for (auto &p : dst) {
          p.x()=recvBuf[i++];
          p.y()=recvBuf[i++];
          p.z()=recvBuf[i++];
          if (format==XYZQ)
            p.charge=recvBuf[i++];
          if (format==XYZQI) {
            p.charge=recvBuf[i++];
            p.id=(typename particle::Tid)recvBuf[i++];
          }
        }
        for (auto &x : recvExtra)
          x=recvBuf[i++];
        assert( (size_t)i==recvBuf.size() );
        if ( (size_t)i!=recvBuf.size() )
          std::cerr << "Particle transmitter says: !!!!!!!!!!!" << endl;
      }

    /*
     * @brief Sum tables computed by parallel processes
     *
     * @details Slave processes send histograms to the master. The master computes the 
     * sum and sends it back to the slaves.
     */
    template<class Ttable>
      void mergeTables(MPIController* mpiPtr, FloatTransmitter &ft, Ttable &table, int &size) {
        if (!mpiPtr->isMaster()) {
          vector<FloatTransmitter::floatp> sendBuf = table.hist2buf(size);
          vector<FloatTransmitter::floatp> recvBuf = ft.swapf(*mpiPtr, sendBuf, mpiPtr->rankMaster());
          table.buf2hist(recvBuf);
        }
        if (mpiPtr->isMaster()) {
          vector<FloatTransmitter::floatp> sendBuf = table.hist2buf(size);
          vector<FloatTransmitter::floatp> recvBuf(size);
          for (int i=0; i<mpiPtr->nproc(); ++i) {
            if (i==mpiPtr->rankMaster()) continue;
            ft.recvf(*mpiPtr, i, recvBuf);
            ft.waitrecv();
            sendBuf.insert(sendBuf.end(), recvBuf.begin(), recvBuf.end());
          }
          table.buf2hist(sendBuf);
          sendBuf = table.hist2buf(size);
          for (int i=0; i<mpiPtr->nproc(); ++i) {
            if (i==mpiPtr->rankMaster()) continue;
            ft.sendf(*mpiPtr, sendBuf, i);
            ft.waitsend();
          }
        }
      }

  } //end of mpi namespace
}//end of faunus namespace

#endif
#endif
