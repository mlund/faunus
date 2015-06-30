#define DIPOLEPARTICLE
#include <faunus/faunus.h>
#include <faunus/multipole.h>
#include <functional>
#include <iostream>
using namespace Faunus;                     
using namespace Faunus::Move;
using namespace Faunus::Potential;

typedef Space<Geometry::Cuboid,DipoleParticle> Tspace; 
typedef DipoleDipole Tpair;

int main() {
	InputMap in("nemo.json");                // open parameter file for user input
	Tspace spc(in);                                // sim.space, particles etc.

	Energy::NonbondedVector<Tspace,Tpair> pot(in); // non-bonded only
	Move::Propagator<Tspace,false> mv(in,pot,spc);

	double resolution = 0.1;
	Histogram<double,unsigned int> HM_x,HM_y,HM_z;
	HM_x.setResolution(resolution);
	HM_y.setResolution(resolution);
	HM_z.setResolution(resolution);
	Table2D<double,Average<double> > M_x, M_y, M_z, En;
	M_x.setResolution(resolution);
	M_y.setResolution(resolution);
	M_z.setResolution(resolution);
	En.setResolution(resolution);

	double sigma = 2.8863;

	spc.p[0] = Point(0,0,6.324*sigma);
	spc.p[0].charge = 0.0;
	spc.p[0].muscalar = 0.0;
	spc.p[0].mu = Point(1.0,0.0,0.0);

	spc.p[2] = Point(0.0,0.0,0.0);
	spc.p[2].charge = 0.0;
	spc.p[2].muscalar = 1.0;
	spc.p[2].mu = Point(1.0,0.0,0.0);
	spc.p[1] = Point(0,0,5.0*sigma);
	spc.p[1].charge = 0.0;
	spc.p[1].muscalar = 1.0;
	spc.p[1].mu = Point(1.0,0.0,0.0);

	for(int i = 3; i < spc.p.size(); i++) {
		spc.p[i] = Point(0.0,0.0,10.0*sigma + double(i));
		spc.p[i].charge = 0.0;
		spc.p[i].muscalar = 0.0;
	}

	spc.trial = spc.p;
	EnergyDrift sys;                               // class for tracking system energy drifts
	sys.init( Energy::systemEnergy(spc,pot,spc.p)  );// initial energy
	
	int cnt = 1;
	MCLoop loop(in);  
	                             // class for mc loop counting
	while ( loop[0] ) {                            // Markov chain 

	double dist = spc.p[1].z();
		while ( loop[1] ) {
			sys += mv.move();

			HM_x(spc.p[2].mu.x()*spc.p[2].muscalar)++;
			HM_y(spc.p[2].mu.y()*spc.p[2].muscalar)++;
			HM_z(spc.p[2].mu.z()*spc.p[2].muscalar)++;
			HM_x(spc.p[1].mu.x()*spc.p[1].muscalar)++;
			HM_y(spc.p[1].mu.y()*spc.p[1].muscalar)++;
			HM_z(spc.p[1].mu.z()*spc.p[1].muscalar)++;
			M_x(dist) += spc.p[2].mu.x()*spc.p[2].muscalar;
			M_y(dist) += spc.p[2].mu.y()*spc.p[2].muscalar;
			M_z(dist) += spc.p[2].mu.z()*spc.p[2].muscalar;
			M_x(dist) += spc.p[1].mu.x()*spc.p[1].muscalar;
			M_y(dist) += spc.p[1].mu.y()*spc.p[1].muscalar;
			M_z(dist) += spc.p[1].mu.z()*spc.p[1].muscalar;
			En(dist) += Energy::systemEnergy(spc,pot,spc.p);
		}    

		spc.p[1] = Point(0,0,5.0*sigma -  ( double(cnt)/10.0 )*4.0*sigma   );
		spc.trial = spc.p;

		cout << "Mu0: " << spc.p[2].mu.transpose() << ", dist: " << dist << ", Mu1: " << spc.p[1].mu.transpose() << endl;
		cnt++;
		sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current
		cout << loop.timing() << std::flush;
	}

	HM_x.save("hist_dip_x.dat");
	HM_y.save("hist_dip_y.dat");
	HM_z.save("hist_dip_z.dat");
	M_x.save("dipdip_x.dat");
	M_y.save("dipdip_y.dat");
	M_z.save("dipdip_z.dat");
	En.save("energy.dat");

	std::cout << spc.info() + pot.info() + mv.info() + sys.info(); // final info

	return 0;
}
