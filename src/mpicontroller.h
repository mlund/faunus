#pragma once
#include "random.h"
#include "core.h"

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cstdio>

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

#define COUT if (isMaster()) std::cout

/**
 * @todo:
 *
 * This is legacy code that works but would benefit from a bit of polishing.
 * In particular:
 * - add exceptions
 * - update documentation
 */

namespace Faunus {

    /**
     * @brief Namespace for Message Parsing Interface (MPI) functionality
     */
    namespace MPI {

        extern std::string prefix;

        /**
         * @brief Main controller for MPI calls
         *
         * This is the MPI controller required for all MPI programs.
         *
         *     MPIController mpi;
         *     mpi.init();
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
                void init();            //!< Initialize MPI and file IO
                void finalize();        //!< Finalize MPI
                ~MPIController();       //!< End of all MPI calls!
                int nproc() const;      //!< Number of processors in communicator
                int rank() const;       //!< Rank of process
                int rankMaster() const; //!< Rank number of the master
                bool isMaster() const;  //!< Test if current process is master
                std::ostream& cout();
                Random random;          //!< Random number generator for MPI calls
                std::string id;         //!< Unique name associated with current rank
#ifdef ENABLE_MPI
                MPI_Comm comm=MPI_COMM_WORLD;    //!< Communicator (Default: MPI_COMM_WORLD)
#endif
            private:
                std::ofstream f; //!< Redirect stdout to here for rank-based file output
                int _nproc=1;      //!< Number of processors in communicator
                int _rank=0;       //!< Rank of process
                int _master=0;     //!< Rank number of the master
                bool mpi_initialized=false;
        };

        void to_json(json&, const MPIController&);

        extern MPIController mpi;

#ifdef ENABLE_MPI

        /**
         * @brief Split N items into nproc parts
         *
         * This returns a pair with the first and last
         * item for the current rank.
         */
        template<class T=int>
            std::pair<T,T> splitEven(const MPIController &mpi, T N) {
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
        double reduceDouble(MPIController &mpi, double local);

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
                std::vector<floatp> swapf(MPIController&, std::vector<floatp>&, int); //!< Swap data with another process
                void sendf(MPIController&, std::vector<floatp>&, int); //!< Send vector of floats
                void recvf(MPIController&, int, std::vector<floatp>&); //!< Receive vector of floats
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
         *     Tpvec myparticles(200); // we have 200 particles
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
        template<typename Tpvec>
            class ParticleTransmitter : public FloatTransmitter {
                private:
                    typedef typename Tpvec::value_type particle; 
                public:
                    enum dataformat {XYZ=3, XYZQ=4, XYZQI=5};
                    std::vector<floatp> sendExtra;                      //!< Put extra data to send here.
                    std::vector<floatp> recvExtra;                      //!< Received extra data will be stored here

                    ParticleTransmitter();
                    void send(MPIController&, const Tpvec&, int); //!< Send particle vector to another node 
                    void recv(MPIController&, int, Tpvec&);       //!< Receive particle vector from another node
                    void waitrecv();
                    void setFormat(dataformat);
                    void setFormat(const std::string&);
                    dataformat getFormat() const;

                private:
                    dataformat format;                             //!< Data format to send/receive - default is XYZQ
                    std::vector<floatp> sendBuf, recvBuf;
                    Tpvec *dstPtr;  //!< pointer to receiving particle vector
                    void pvec2buf(const Tpvec&); //!< Copy source particle vector to send buffer
                    void buf2pvec(Tpvec&);       //!< Copy receive buffer to target particle vector
            };

        template<typename Tpvec>
            ParticleTransmitter<Tpvec>::ParticleTransmitter() { setFormat(XYZQI); }

        template<typename Tpvec>
            void ParticleTransmitter<Tpvec>::setFormat(dataformat d) { format = d; }

        template<typename Tpvec>
            void ParticleTransmitter<Tpvec>::setFormat(const std::string &s) {
                setFormat(XYZQI);
                if (s=="XYZQ")
                    setFormat(XYZQ);
                if (s=="XYZ")
                    setFormat(XYZ);
            }

        template<typename Tpvec>
            typename ParticleTransmitter<Tpvec>::dataformat
            ParticleTransmitter<Tpvec>::getFormat() const { return format; }

        /*!
         * \param mpi MPI controller to use
         * \param src Source particle vector
         * \param dst Destination node
         */
        template<typename Tpvec>
            void ParticleTransmitter<Tpvec>::send(MPIController &mpi, const Tpvec &src, int dst) {
                assert(dst>=0 && dst<mpi.nproc() && "Invalid MPI destination");
                pvec2buf(src);
                FloatTransmitter::sendf(mpi, sendBuf, dst);
            }

        template<typename Tpvec>
            void ParticleTransmitter<Tpvec>::pvec2buf(const Tpvec &src) {
                sendBuf.clear();
                for (auto &p : src) {
                    sendBuf.push_back(p.pos.x());
                    sendBuf.push_back(p.pos.y());
                    sendBuf.push_back(p.pos.z());
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
        template<typename Tpvec>
            void ParticleTransmitter<Tpvec>::recv(MPIController &mpi, int src, Tpvec &dst) {
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

        template<typename Tpvec>
            void ParticleTransmitter<Tpvec>::waitrecv() {
                FloatTransmitter::waitrecv();
                buf2pvec(*dstPtr);
            }

        template<typename Tpvec>
            void ParticleTransmitter<Tpvec>::buf2pvec(Tpvec &dst) {
                int i=0;
                for (auto &p : dst) {
                    p.pos.x()=recvBuf[i++];
                    p.pos.y()=recvBuf[i++];
                    p.pos.z()=recvBuf[i++];
                    if (format==XYZQ)
                        p.charge=recvBuf[i++];
                    if (format==XYZQI) {
                        p.charge=recvBuf[i++];
                        p.id=(int)recvBuf[i++];
                    }
                }
                for (auto &x : recvExtra)
                    x=recvBuf[i++];
                assert( (size_t)i==recvBuf.size() );
                if ( (size_t)i!=recvBuf.size() )
                    std::cerr << "Particle transmitter says: !!!!!!!!!!!" << std::endl;
            }

        /*
         * @brief Sum tables computed by parallel processes
         *
         * @details Slave processes send histograms to the master. The master computes the 
         * average and sends it back to the slaves. Ttable can be Table, Table2D or Table3D in auxiliary.h.
         *
         */
        template<class Ttable>
            void avgTables(MPIController* mpiPtr, FloatTransmitter &ft, Ttable &table, int &size) {
                if (!mpiPtr->isMaster()) {
                    std::vector<FloatTransmitter::floatp> sendBuf = table.hist2buf(size);
                    std::vector<FloatTransmitter::floatp> recvBuf = ft.swapf(*mpiPtr, sendBuf, mpiPtr->rankMaster());
                    table.buf2hist(recvBuf);
                }
                if (mpiPtr->isMaster()) {
                    std::vector<FloatTransmitter::floatp> sendBuf = table.hist2buf(size);
                    std::vector<FloatTransmitter::floatp> recvBuf(size);
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
#endif

    } //end of mpi namespace
}//end of faunus namespace

