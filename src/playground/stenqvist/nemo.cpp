#include <faunus/faunus.h>

// Manual: http://faunus.sourceforge.net/doxyhtml/index.html

using namespace Faunus;
using namespace std;

class Point {
  public:
    double x,y,z;
    void init(double,double,double);
};

class particle : public Point {
  public:
    std::string name;
    double sigma;
    double eps;
    double q;
    Point mu;
    Eigen::Matrix3d alpha;
    void init(double,double,double,double,double,double,Point,Eigen::Matrix3d);
};

void particle::init (double sigma_in, double eps_in, double x_in, double y_in, double z_in, double q_in, Point mu_in,Eigen::Matrix3d alpha_in) {
  sigma = sigma_in;
  eps = eps_in;
  x = x_in;
  y = y_in;
  z = z_in;
  q = q_in*pc::charge;
  mu = mu_in;
  alpha = alpha_in;
}

void Point::init (double x_in, double y_in, double z_in) {
  x = x_in;
  y = y_in;
  z = z_in;
}

double dist(const Point &p1, const Point &p2) {
  double dx = p1.x - p2.x;
  double dy = p1.y - p2.y;
  double dz = p1.z - p2.z;
  return sqrt(dx*dx+dy*dy+dz*dz);
}

Point getRab(const Point &p1, const Point &p2) {
  Point p;
  p.x = p2.x - p1.x;
  p.y = p2.y - p1.y;
  p.z = p2.z - p1.z;
  return p;
}

double DP(const Point &p1, const Point &p2) {
  return (p1.x*p2.x + p1.y*p2.y + p1.z*p2.z);
}

Point getUnitVector(const particle p1, const particle p2) {
  Point R = getRab(p1,p2);
  double r = dist(p1,p2);
  R.x = R.x/r;
  R.y = R.y/r;
  R.z = R.z/r;
  return R;
}

Point rotatePoint(Point p, double dTheta, double dPhi) {
  double x,y,z;
  Point q = getSPC(p);
  q.y = q.y + dTheta;
  q.z = q.z + dPhi;
  x = q.x*std::sin(q.y)*std::cos(q.z);
  y = q.x*std::sin(q.y)*std::sin(q.z);
  z = q.x*std::cos(q.y);
  p.x = x;
  p.y = y;
  p.z = z;
  return p;
}

double calcKeesom(const particle p1, const particle p2) {
  //std::cout << " Correct? " << DP(p1.mu,p2.mu)/(4*pc::den*pow(dist(p1,p2),3)) << " < " << pc::kT << "\n";
  return -DP(p1.mu,p1.mu)*DP(p2.mu,p2.mu)/(3*pc::den*pc::den*pc::kT*pow(dist(p1, p2),6));
}

double calcDIMean(const particle p1, const particle p2) {
  return -p1.q*p1.q*DP(p2.mu,p2.mu)/(6*pc::den*pc::den*pc::kT*pow(dist(p1,p2),4));
}

double getRand() {
  return ((double)rand()/RAND_MAX);
}

double getRandn() {
  return (getRand()-0.5)/0.5;
}

Point ranunit() {
  Point u;
  double r2;
  do {
    u.x=2*( getRand()-0.5 );
    u.y=2*( getRand()-0.5 );
    u.z=2*( getRand()-0.5 );
    r2=u.x*u.x+u.y*u.y+u.z*u.z;
  } while (r2>1);
  r2=sqrt(r2);
  u.x=u.x/r2;
  u.y=u.y/r2;
  u.z=u.z/r2;
  return u;
}

Point randomRotation(Point p, double step1, double step2) {
  Point u=ranunit();
  Eigen::Vector3d a;
  double theta=getRandn()*pc::pi/2;
  a.x()=u.x;
  a.y()=u.y;
  a.z()=u.z;
  Eigen::Quaterniond q;
  q=Eigen::AngleAxisd(theta, a);
  a.x()=p.x;
  a.y()=p.y;
  a.z()=p.z;
  a=q*a;
  p.x=a.x();
  p.y=a.y();
  p.z=a.z();
  return p;
}

std::vector<particle> getStartConfiguration(double r) {
  particle par1,par2;
  Point mu;
  std::vector<particle> p(2);

  Eigen::Matrix3d alpha(3,3); // se "initializer lists"
  alpha(0,0) = 7.14264826;
  alpha(1,0) = 0.00014943;
  alpha(1,1) = 5.64937691;
  alpha(2,0) = -0.00003816;
  alpha(2,1) = 0.00035990;
  alpha(2,2) = 6.03919669;

  mu.init(0.0,0.0,2.4126*pc::debye);
  par1.init(0.8, 0.1, 0, 0, 0, -1.0, mu,alpha);
  p.at(0) = par1;

  mu = randomRotation(mu,100,100);
  par2.init(0.8, 0.1, r, 0, 0, +1.0, mu,alpha);
  p.at(1) = par2;

  return p;
}

particle randomRadius(particle p, double dr) {
  p.x += dr*getRandn();
  return p;
}

std::vector<particle> initialRandom(std::vector<particle> p, double step, double dr, int opt) {
  double a = getRand();

  switch(opt) {
    case 1:
      if(a < 0.5) {
        p.at(0).mu = randomRotation(p.at(0).mu,step,step);
      } else {
        p.at(1).mu = randomRotation(p.at(1).mu,step,step);
      }
      break;
    case 2:
      p.at(1).mu = randomRotation(p.at(1).mu,step,step);
      break;
    case 3:
      if(a < 0.333333) {
        p.at(0).mu = randomRotation(p.at(0).mu,step,step);
      } else if(a < 0.666667) {
        p.at(1).mu = randomRotation(p.at(1).mu,step,step);
      } else {
        p.at(1) = randomRadius(p.at(1),dr);
      }
      break;
    case 4:
      if(a < 0.5) {
        p.at(1) = randomRadius(p.at(1),dr);
      } else {
        p.at(1).mu = randomRotation(p.at(1).mu,step,step);
      }
      break;
    case 5:
      p.at(1) = randomRadius(p.at(1),dr);
  }
  return p;
}

std::vector<particle> checkEnergy(double dE,std::vector<particle> p,std::vector<particle> p_t) {
  /*double r = dist(p_t.at(0),p_t.at(1));
    if (r < 2.0e-10 || r > 20e-10) {
    return p;
    }*/

  if(getRand() > std::exp(-dE) ) {
    return p;
  }
  accepted++;
  return p_t;
}

double calculateEnergy(std::vector<particle> &p, Energy &ene, int opt) {
  double r = dist(p.at(0),p.at(1));
  if (r < 0.5e-10 || r > 10e-10) {
    return 1e9;
  }

  if (opt == 2 || opt == 4) {
    return ene.q2mu(p.at(0), p.at(1))/(pc::kT);
  } else if (opt == 1 || opt == 3) {
    return ene.mu2mu(p.at(0), p.at(1))/(pc::kT);
  } else if (opt == 5) {
    return ene.q2q(p.at(0),p.at(1))/(pc::kT);
  } else {
    return 0.0;
  }
}

double getMean(std::vector<particle> p, int opt) {
  if(opt == 1 || opt == 3) {
    return calcKeesom(p.at(0), p.at(1))/pc::kT;
  } else if (opt == 2 || opt == 4) {
    return calcDIMean(p.at(0), p.at(1))/pc::kT;
  } else {
    return 0.0;
  }
}

class DipoleParticle : public PointParticle {
  public:
    double mu_s;            //!< Dipole moment scalar
    Eigen::Vector3d mu;     //!< Dipole moment unit vector
    Eigen::Matrix3d alpha;  //!< Polarization tensor

    inline DipoleParticle() {
      mu_s=0;
    }

    template<typename OtherDerived>                                   
      DipoleParticle(const Eigen::MatrixBase<OtherDerived>& other) : Tvec(other) {}  

    template<typename OtherDerived>                                   
      DipoleParticle& operator=(const Eigen::MatrixBase<OtherDerived> &other) {      
        Tvec::operator=(other);                                       
        return *this;                                                 
      }  
};

Point getSPC(Point p) {
  double r,theta,phi;
  r = std::sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
  if (r != 0) {
    theta = std::acos(p.z/r);
    if (p.x != 0) {
      phi = std::atan(p.y/p.x);
    } else {
      if (std::sin(theta) != 0) {
        phi = std::acos(r*std::sin(theta));
      } else {
        phi = 0;
      }
    }
  } else {
    theta = 0;
    phi = 0;
  }
  Point q;
  q.x = r;
  q.y = theta;
  q.z = phi;
  return q;
}

class Energy {
  public:
    double q2q(const particle &p1, const particle &p2);
    double q2mu(const particle &p1, const particle &p2);
    double q2Q(const particle &p1, const particle &p2);
    double q2alpha(const particle &p1, const particle &p2); 
    double mu2mu(const particle &p1, const particle &p2); 
    double mu2Q(const particle &p1, const particle &p2); 
    double mu2alpha(const particle &p1, const particle &p2); 
    double Q2Q(const particle &p1, const particle &p2);
    double Q2alpha(const particle &p1, const particle &p2); 
    double alpha2alpha(const particle &p1, const particle &p2); 
    double total(const particle &p1, const particle &p2);
};

double Energy::q2q(const particle &p1, const particle &p2) {
  double R = dist(p1,p2);
  return p1.q*p2.q/(R*pc::den);
}

double Energy::q2mu(const particle &p1, const particle &p2) {
  Point r_0;
  r_0.init(0.0,0.0,0.0);
  Point R = geo.vdist(p1,p2);
  //Point R = getRab(p1,p2);
  double cosT = DP(R,p2.mu)/(dist(R,r_0)*dist(p2.mu,r_0));
  Point q = getSPC(p2.mu);
  double r = dist(p1,p2);
  return -p1.q*q.x*cosT/(pc::den*r*r);
}

double Energy::mu2mu(const particle &p1, const particle &p2) {
  double R = dist(p1,p2);
  Point q1 = getSPC(p1.mu);
  Point q2 = getSPC(p2.mu);
  return -q1.x*q2.x*(2*std::cos(q1.y)*std::cos(q2.y)-std::sin(q1.y)*std::sin(q2.y)*std::cos(q1.z-q2.z))/(pc::den*pow(R,3));
}

// Insert energy functions here!
// Major goal:
//   function: pairpot(particle &a, particle &b, point &vdist)

int main() {
  InputMap mcp("nemo.conf");
  Geometry::Sphere geo(mcp);
  Potential::Coulomb coulomb(mcp);

  cout << geo.info() << coulomb.info(25);

  DipoleParticle a, b;

  Point vdist = geo.vdist(a,b);
  double u = coulomb(a,b, geo.dist(a,b) );

  a.ranunit( slp_global ); // random unit vector
  double r = slp_global(); // random number [0:1[

  // Examples
  cout << "Squared distance = " << geo.sqdist(a,b) << "\n"; 
  cout << "Bjerrum          = " << coulomb.bjerrumLength() << "\n";
  cout << "Avogadro's number= " << pc::Nav << "\n";  // se "faunus/physconst.h"
}

int main1() {
  accepted=0;
  double r, dr,limit,step,dE;		// Initialize variables
  std::vector<particle> p,p_t;
  Energy ene;
  int N = 1e7;					// Number of random displacements
  r = 10e-10;					    // Radius between particles
  dr = 9e-10;						// Size of radius displacement
  step = 100.0;					// % of possible angle displacement to be displaced
  std::vector<Point> mu1(1);
  std::vector<Point> mu2(1);

  double dbin=0.1;
  std::vector<int> hist(1000);

  for (int i=0; i < hist.size(); i++) {
    hist[i]=0;
  }

  int opt = 3;					// 1 = dipole-dipole (no dr change)      2 = ion-dipole (no dr change)
  // 3 = dipole-dipole (include dr change) 4 = ion-dipole (include dr change)
  // 5 = charge-charge (include dr change)

  p = getStartConfiguration(r);										// Get initial configuration

  for(int n = 1; n < N ; n++) {
    p_t = initialRandom(p,step,dr,opt);								// Random displacement

    dE = calculateEnergy(p_t,ene,opt) - calculateEnergy(p,ene,opt); // Calculate energy differenace
    p = checkEnergy(dE,p,p_t);				                        // Check if move is allowed

    double r=dist(p[0],p[1])*1e10;
    int bin=(int)( (r/dbin)+0.5 );
    hist[bin]++;
  }

  std::ofstream oh("tabHist.dat");
  for (int i=0; i<hist.size(); i++)
    if (hist[i]>0)
      oh << i*dbin << " " << hist[i] << "\n";

  cout << "Acceptance = " << accepted/double(N) << "\n";

  //endLine();
  return 1;
}
