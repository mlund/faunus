/* 
You are only allow to use this software if you have signed the following license. 

ACADEMIC/NON-PROFIT 
SASA SOFTWARE LICENSE AGREEMENT
IMPORTANT: This SASA license Agreement is a legal agreement between you, the end user (either an individual or an entity), and the Karlsruhe Institute of Technology. 
SASA Software License
GRANT OF LICENSE. Karlsruhe Institute of Technology grants, and you hereby accept, a non-exclusive license to use the SASA software product of the version specified above ("Software") to the extent of its rights and in accordance with the terms of this Agreement. This licensed copy of the Software may only be used on computers at your site by you and members of your organization at your site who have read and agreed to this license. You may install the Software on computers at your site for your own use or use by members of your organization at your site. You may not distribute copies of the Software to others outside of your site. You may make only those copies of the Software which are necessary to install and use it as permitted by this Agreement, or are for purposes of backup and archival records. 
OWNERSHIP. This ownership is protected by the copyright laws of the Federal Republic of Germany and by international treaty provisions. Upon expiration or termination of this Agreement, you shall promptly return all copies of the Software and accompanying written materials to the Karlsruhe Institute of Technology. 
MODIFICATIONS AND DERIVATIVE WORKS. You may modify the software, and use it to create derivative works, for your internal use at the site covered by this license. You may not distribute such modified or derivative software to others outside of your site without written permission. You may distribute the modifications themselves (e.g. as "patches") under terms of your choice. We encourage users to contribute modifications back into the Software, but you are under no obligation to do so. 
REPORTS OF PUBLICATIONS. You agree to acknowledge use of the Software in any reports or publications of results obtained with the Software and cite the publications listed on the download page where you obtained the software in any report or publication in which the Software was used. If you fail to properly acknowledge the use of the software you agree to pay the industrial license fee.
ASSIGNMENT RESTRICTIONS. You shall not use the Software (or any part thereof) in connection with the provision of consultancy, modeling or other services, whether for value or otherwise, on behalf of any third party who does not hold a current valid SASA  Software License Agreement. You shall not use the Software to write other software that duplicates the functionality of the Software. You shall not rent, lease, or otherwise sublet the Software or any part thereof. You may transfer on a permanent basis the rights granted under this license provided you transfer this Agreement and all copies of the Software, including prior versions, and all accompanying written materials. The recipient must agree to the terms of this Agreement in full and register this transfer with the Karlsruhe Institute of Technology. 
LIMITED WARRANTY. LICENSEE acknowledges that LICENSORS make no warranty, expressed or implied, that the program will function without error, or in any particular hardware environment, or so as to generate any particular function or result, and further excluding any other warranty, as to the condition of the program, its merchantability, or its fitness for a particular purpose. LICENSORS shall not be liable for any direct, consequential, or other damages suffered by the LICENSEE or any others as a result of their use of the program, whether or not the same could have been foreseen by LICENSORS prior to granting this License. In no event shall LICENSORS liability for any breach of this agreement exceed the fee paid for the license. 
KARLSRUHE INSTITUTE OF TECHNOLOGY'S LIABILITY. In no event shall the Karlsruhe Institute of Technology be liable for any indirect, special, or consequential damages, such as, but not limited to, loss of anticipated profits or other economic loss in connection with or arising out of the use of the software by you or the services provided for in this Agreement, even if the Karlsruhe Institute of Technology has been advised of the possibility of such damages. The Karlsruhe Institute of Technology's entire liability and your exclusive remedy shall be, at the Karlsruhe Institute of Technology's discretion, to return the Software and proof of purchase to the Karlsruhe Institute of Technology for either (a) return of any license fee, or (b) correction or replacement of Software that does not meet the terms of this limited warranty. 
NO OTHER WARRANTIES. The Karlsruhe Institute of Technology disclaims other implied warranties, including, but not limited to, implied warranties of merchantability or fitness for any purpose, and implied warranties arising by usage of trade, course of dealing, or course of performance. Some states do not allow the limitation of the duration or liability of implied warranties, so the above restrictions might not apply to you. 
LICENSE FEE. The software is free for non-profit, government and academic organizations. For-profit and commercial organizations wishing to license SASA shall contact: 
SASA-support@kit.edu 
and request a quote for a commercial license.
Supercomputer centers can license SASA under the same conditions and make it available to their users from non-profit organizations as executable code. However, for-profit organizations who want to use the program at supercomputer centers must have signed a separate license agreement.


If you have no license please contact SASA-support@kit.edu
*/


#ifndef POWERSASA_H_
#define POWERSASA_H_

// Uncomment the next line to get volumes of individual atoms
#define ATOMVOL

#include <cstdlib>
#include <cmath>
#include <iostream>
#include "power_diagram.h"
#include <vector>

#define MAX_NB 20
#define MAX_VX 12
#define MAX_PNT 4
#define MAX_COUNT 100
//=============================================================

namespace POWERSASA
{


class PowerSasaException : public std::exception {};

/*!
 * \class PowerSasa
 * \brief Sasa calculation class
 *
 * This class calculates the Sasa and stores it by using a power diagram.
 *
 * \endcode
 */

template <class PDFloat, class PDCoord>      // floating type and vector type
class PowerSasa 
{
public:
	inline static PDFloat DRAD2() { return static_cast<PDFloat>(1000.0) * std::numeric_limits<PDFloat>::epsilon() ; }
	inline static PDFloat DANG() { return static_cast<PDFloat>(1000.0) * std::numeric_limits<PDFloat>::epsilon() ; }
	inline static PDFloat pi() { return static_cast<PDFloat> (3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342) ;} //everything const, will be optimized away.

	const std::vector<PDFloat>& getSasa() const {return Sasa;}
	const std::vector<PDFloat>& getVol() const {return Vol;}
	const std::vector<PDCoord>& getDVol() const {return DVol;}
	const std::vector<PDCoord>& getDSasa() const {return DSasa;}

	void calc_sasa_single(const unsigned int iatom);
	void calc_sasa_all();

	template<class Coordcontainer, class Floatcontainer>
	void update_coords(Coordcontainer const& coords, Floatcontainer const& radii)
	{
		unsigned int n_old = power_diagram->get_points().size();
		power_diagram->recalculate(coords.begin(),radii.begin(),coords.size());
		if (n_old < power_diagram->get_points().size()) Resize_NA();
	}

	template <class Pos_iterator, class Strength_iterator>
	void add_more(const Pos_iterator pos_it, const Strength_iterator strength_it, const unsigned int newSize) {
		power_diagram->addMore(pos_it, strength_it, newSize);
		Resize_NA();
	}

	void add_more(const PDCoord& position, const PDFloat& radius,const int near) {
		power_diagram->addMore(position,radius,near);
		Resize_NA();
	}
	void revert() {power_diagram->revert();}
	POWER_DIAGRAM::PowerDiagram<PDFloat,PDCoord,3 > & get_power_diagram()  { return *power_diagram; }

	unsigned int NumOfNeighbours(unsigned int iatom) const
	{
		return power_diagram->get_points()[iatom].neighbours.size();
	}
	unsigned int AtomNo(unsigned int i_atom, unsigned int i_neighbour) const
	{
		return (unsigned int)(power_diagram->get_points()[i_atom].neighbours[i_neighbour]
			- &power_diagram->get_points()[0]);
	}

	template<class Coordcontainer, class Floatcontainer, class Intcontainer>
	PowerSasa(Coordcontainer const& coords, Floatcontainer const& radii, Intcontainer const& bond_to,
		const bool with_Sasa, const bool with_dSasa, const bool with_Vol, const bool with_dVol) :
		withSasa(with_Sasa), withDSasa(with_dSasa), withVol(with_Vol), withDVol(with_dVol), power_diagram(0)
	{
		power_diagram = new POWER_DIAGRAM::PowerDiagram<PDFloat, PDCoord,3>(POWER_DIAGRAM::PowerDiagram<PDFloat, PDCoord,3>::create(coords.size(),coords.begin(),radii.begin(),bond_to.begin())
			.with_radiiGiven(1).with_calculate(1).with_cells(1).with_zeroPoints(1).with_Warnings(0).withoutCheck(1));
		Init();
	}
	template<class Coordcontainer, class Floatcontainer>
	PowerSasa(Coordcontainer const& coords, Floatcontainer const& radii,
			const bool with_Sasa=1, const bool with_dSasa=0, const bool with_Vol=0, const bool with_dVol=0) :
			power_diagram(0), withSasa(with_Sasa), withDSasa(with_dSasa), withVol(with_Vol), withDVol(with_dVol)
	{
		std::vector<int> bond_to;
		bond_to.reserve(coords.size());
		bond_to.push_back(0);
		for (unsigned int i = 1; i < coords.size(); ++i)
		{
			bond_to.push_back(i-1);
		}
		power_diagram = new POWER_DIAGRAM::PowerDiagram<PDFloat, PDCoord,3>(POWER_DIAGRAM::PowerDiagram<PDFloat, PDCoord,3>::create(coords.size(),coords.begin(),radii.begin(),bond_to.begin())
		        .with_radiiGiven(1).with_calculate(1).with_cells(1).with_zeroPoints(1).with_Warnings(0).without_Check(1));
		Init();
	}

	virtual ~PowerSasa() 
	{
		if (power_diagram != 0)
		{
			delete power_diagram;
			power_diagram = 0;
		}
	}
private:
	inline void Init()
	{
		Resize_NA();
		Resize_NB(MAX_NB);
		Resize_VX(MAX_VX);
		Resize_PNT(MAX_PNT);
	}
	inline void Resize_NB(unsigned int nnb)
	{
		np.resize(nnb);
		nt.resize(nnb);
		e.resize(nnb);
		sintheta.resize(nnb);
		costheta.resize(nnb);
		nb_RAD2.resize(nnb);
		nb_dist.resize(nnb);
#ifdef ATOMVOL
		volnb.resize(nnb);
		knot.resize(nnb);
		fknot.resize(nnb);
#endif
	      	unsigned int npnt = MAX_PNT;
		if (p.size() > 0) npnt = p[0].size();
		unsigned int nnb_old = p.size();
		next.resize(nnb);
		p.resize(nnb);
		ang.resize(nnb);
		for (unsigned int i = nnb_old; i < nnb; ++i)
		{
			next[i].resize(npnt);
			p[i].resize(npnt);
			ang[i].resize(npnt);
		}
		if (withDSasa)
		{
		        for (unsigned int i = 0; i < DSasa_parts.size(); ++i) DSasa_parts[i].resize(nnb);
		}
	}
        inline void Resize_VX(unsigned int nvx)
	{
		off.resize(nvx);
		vx.resize(nvx);
		br_c.resize(nvx);
		br_p.resize(nvx);
		for (unsigned int i = 0; i < nvx; ++i)
		{
			br_c[i].resize(2);
			br_p[i].resize(2);
		}
	}
	inline void Resize_PNT(unsigned int npnt)
	{
		for (unsigned int i = 0; i < p.size(); ++i)
		{
			next[i].resize(npnt);
			p[i].resize(npnt);
			ang[i].resize(npnt);
		}
	        rang.resize(npnt);
		pos.resize(npnt);
	}
	inline void Resize_NA()
	{
		unsigned int n = power_diagram->get_points().size();
		if (withSasa) Sasa.resize(n,0);
		if(withDSasa)
		{
			unsigned int nnb = MAX_NB;
		        if(DSasa_parts.size()>0)if (DSasa_parts.front().size() > 0) nnb = DSasa_parts.front().size();
		        unsigned int n_old = DSasa.size();
			DSasa.resize(n);
			DSasa_parts.resize(n);
		        for (unsigned int i = n_old; i < n; ++i) DSasa_parts[i].resize(nnb);
			
		}
		if (withVol) Vol.resize(n,0);
		if (withDVol) DVol.resize(n,PDCoord(0,0,0));

		PDFloat maxr2 = 0.0;
		for (unsigned int i = 0; i < n; ++i)
		{
			if (power_diagram->get_points()[i].r2 > maxr2)
			  maxr2 = power_diagram->get_points()[i].r2;
		}
		tol_pow = maxr2 * DRAD2();
	}

 	inline PDFloat Ang_About(PDCoord const& a, PDCoord const& b, PDCoord const& c);
	inline void Get_Ang(const int & np, const std::vector<int> & p, const PDCoord & e,
			const PDFloat & sintheta, const PDFloat & costheta, std::vector<PDFloat> & ang);
        inline void Get_Next(int n, std::vector<PDFloat> & ang, std::vector<int> & next,
			const std::vector<int> & p, const PDCoord & e);
	POWER_DIAGRAM::PowerDiagram<PDFloat,PDCoord,3> *power_diagram;
	
	const bool withSasa;
	const bool withDSasa;
	const bool withVol;
	const bool withDVol;
	std::vector<PDFloat>                  Sasa;
	std::vector< std::vector<PDCoord> >   DSasa_parts;
	std::vector<PDCoord>		      DSasa;
	std::vector<PDFloat>                  Vol;
	std::vector<PDCoord>                  DVol;
	
	PDFloat                               tol_pow;

	// ----- for calc_sasa_single --------------------

	std::vector<int> np;            // number of points (registered vertices) of i-th atom
	std::vector<int> nt;            // counts situations that exclude "single cirle"
	std::vector<PDCoord> e;         // direction to neighbour
        std::vector<PDFloat> sintheta;
        std::vector<PDFloat> costheta;
	std::vector<PDFloat> nb_RAD2;
	std::vector<PDFloat> nb_dist;

	std::vector<int> off;           // off[i] != 0 if vertix i is already taken into account
	std::vector<PDCoord> vx;        // all surface vertices
	std::vector< std::vector<int> > br_c;       // bridge between circles (intersections with neighbors)
	std::vector< std::vector<int> > br_p;       // bridge between point numbers
  
	std::vector< std::vector<PDFloat> > ang; // angles
	std::vector< std::vector<int> >next;     // next[][i] is the number of angle that follows the i-th one
	std::vector< std::vector<int> > p;       // number of surface vertices

#ifdef ATOMVOL
	std::vector<PDFloat> volnb;
	std::vector<PDCoord> knot;
	std::vector<bool>    fknot;
#endif
        //----- for Get_Next ----------------------------

	std::vector<int> rang;
	std::vector<int> pos;
	// special thing
public:
	template <const int n>
	void calc_sasa_all(const PDFloat steps[n],PDFloat resultS[][n+1],PDFloat resultV[][n+1])
	{
		for (unsigned int i = 0; i < power_diagram->getPoints().size(); ++i)
		{
			calc_sasa_single(i);
			resultS[i][0]=Sasa[i];
			resultV[i][0]=Vol[i];
		}
		const int& size=power_diagram->getPoints().size();
		for(int a=0;a<size;a++)
		{
			for(int s=0;s<n;s++)
			{
				add_more(power_diagram->getPoints()[a].position+power_diagram->center,power_diagram->getPoints()[a].r+steps[s],a);
				calc_sasa_single(size);
				resultS[a][s]=Sasa[size];
				resultV[a][s]=Vol[size];
				revert();
			}
		}
	};
};

//Implementation:

//-----------------------------------------------------------------------------
// Retuns twist angle (in interval [0, 2*pi()) ) between Vectors a and b about c.
// a and b should be of unit length and perpendicular to c.

template <class PDFloat, class PDCoord> PDFloat PowerSasa<PDFloat, PDCoord>::
Ang_About(PDCoord const& a, PDCoord const& b, PDCoord const& c)
{
	PDFloat const ONE = 0.999;
	PDFloat const THRESHOLD = 0.001;
	PDFloat ang, co, vp;
	PDCoord  v;

	co = a.dot(b);
	if (co <= -ONE)
	{
		v = a.cross(b);
		ang = pi() - std::asin( v.norm() );
	}
	else if (co >= ONE)
	{
		v = a.cross(b);
		ang = std::asin( v.norm() );
	}
	else ang = std::acos(co);

	if (std::fabs(c[0]) > THRESHOLD)
	{
		vp = a[1]*b[2] - a[2]*b[1];
		if ( (vp < 0.0) != (c[0] < 0.0) ) ang = - ang;
	}
	else if (std::fabs(c[1]) > THRESHOLD)
	{
		vp = a[2]*b[0] - a[0]*b[2];
		if ( (vp < 0.0) != (c[1] < 0.0) ) ang = - ang;
	}
	else if (std::fabs(c[2]) > THRESHOLD)
	{
		vp = a[0]*b[1] - a[1]*b[0];
		if ( (vp < 0.0) != (c[2] < 0.0) ) ang = - ang;
	}
	else
	{
		std::cerr << "PowerSasa: Axis too short" << std::endl;
		throw PowerSasaException();
	}
	if (ang < 0) ang += 2.0 * pi();
	return ang;
}
//-----------------------------------------------------------------------------
//Finds phi-angles of s-vertices

template <class PDFloat, class PDCoord> void PowerSasa<PDFloat, PDCoord>::
Get_Ang(const int & np,                              // total number of s-vertices for given circle
        const std::vector<int> & p,                  // s-vertice number for given circle
	const PDCoord & e,                           // direction to neighbour
        const PDFloat & sintheta, const PDFloat & costheta,  // theta angle
	std::vector<PDFloat> & ang)                          // output
	
{
//	if (np == 0) return;//is tested outside
	static PDCoord pu0, pu;
	ang[0] = 0.0;
	pu0 = (vx[p[0]] - costheta * e) / sintheta;
	for (int j = 1; j < np; j++)
	{
		pu = (vx[p[j]] - costheta * e) / sintheta;
		ang[j] =  Ang_About (pu0, pu, e);
	}
}


//-----------------------------------------------------------------------------
// given array ang of length n founds number next[i] of ang element that is next in magnitude to element i
// can be optimized (???)

template <class PDFloat, class PDCoord> void PowerSasa<PDFloat, PDCoord>::
Get_Next(int n, std::vector<PDFloat> & ang, std::vector<int> & next,
	 const std::vector<int> & p, const PDCoord & e)

{
	int j, k, m;
	//static int rang[MAX_PNT], pos[MAX_PNT];
	if (n == 0) return;

	for (j = 1; j < n; ++j)
	{
		if (ang[j] <= 2.0*pi() - DANG()) continue;
		PDFloat dp = vx[p[j]].cross(vx[p[0]]).dot(e);
		if (dp < 0.0) ang[j] = 0.0;
		else if (dp == 0.0)
		{
			std::cout << "PowerSasa: Precision insufficient to resolve angles" << std::endl;
			throw PowerSasaException();
		}
	}

	rang[0] = 0;
	rang[1] = 1;

	for (j = 2; j < n; j++)
	{
		m = j;
		for (k = 1; k < j; k++)
		{
			if (ang[k] > ang[j] + DANG())
			{
				(rang[k])++;
				m--;
			}
			else if (ang[k] > ang[j] - DANG())
			{

			        PDFloat dp = vx[p[j]].cross(vx[p[k]]).dot(e);
				if (dp > 0.0)
				{
					(rang[k])++;
					m--;
			        }
			        else if (dp == 0.0)
				{
					std::cout << "PowerSasa: Precision insufficient to resolve angles" << std::endl;
					throw PowerSasaException();
				}
			}
		}
		rang[j] = m;
	}

	for (j = 0; j < n; j++) pos[rang[j]] = j;

	for (j = 0; j < n; j++)
	{
		if (rang[j] == n-1) next[j] = pos[0];
		else next[j] = pos[ rang[j]+1 ];
	}
}
//===================================================================================

template <class PDFloat, class PDCoord> void PowerSasa<PDFloat, PDCoord>::
calc_sasa_single(const unsigned int iatom)
{

// ----- some initialisations  --------------------

	int i, j, kn, ok;

	std::vector <typename POWER_DIAGRAM::PowerDiagram<PDFloat,PDCoord,3>::cell >
		const &atoms = power_diagram->get_points();
	const typename POWER_DIAGRAM::PowerDiagram<PDFloat,PDCoord,3>::cell	 &atom = atoms[iatom];
	std::vector<typename POWER_DIAGRAM::PowerDiagram<PDFloat,PDCoord,3>::cellPtr> 
		const &nbs = atom.neighbours;
	const int nnb = nbs.size();               // number of neighbors
	PDCoord const &pos = atom.position;       // my coordinates

	if (nnb >= static_cast<int>(np.size()))
	{
		Resize_NB(nnb);
		//std::cerr << "Power Sasa: number of neighbors exeeds MAX_NB=" << MAX_NB << std::endl;
		//throw PowerSasaException();
	}

	if (withSasa != 0) Sasa[iatom] = 0.0;
	if (withDSasa != 0){DSasa[iatom].setZero(); for (i = 0; i < nnb; ++i) DSasa_parts[iatom][i].setZero();}
	if (withVol != 0) Vol[iatom] = 0.0;
	if (withDVol != 0) DVol[iatom].setZero();
	int do_sasa = (withSasa != 0 || withVol != 0) ? 1 : 0;

	PDFloat RAD  = atom.r;
	PDFloat RAD2 = atom.r2;

	if (nnb == 0)
	{
 		if(power_diagram->get_vertices()[0].generators[0]==&atom)
		{
			if (withSasa)Sasa[iatom] = 4*3.1415926535897932384626433832795*RAD2;
			if(withVol)Vol[iatom] = 4*3.1415926535897932384626433832795*0.33333333333333333333333333333333*RAD*RAD2;
		}
		return;
	}

// check that the Power values are not close to zero

	int covered = 1;
	ok = 1;
	
	for (unsigned int n = 0; n < atom.myVertices.size(); ++n)
	{
		if (std::fabs(atom.myVertices[n]->powerValue) < tol_pow)
		{
			std::cerr << "PowerSasa: Vertex power value is too close to zero" << std::endl;
			throw PowerSasaException();
		}
		if (atom.myVertices[n]->powerValue > 0.0) covered = 0;
	}

#ifdef ATOMVOL
	if (covered == 1 && withVol == 0) return;
#else
	if (covered == 1) return;
#endif

//----- get angles and other properties of neighbors -----

	int n_apart = 0;                // number of non-contributing neighbours
	PDFloat dist, nb_RAD;           // distance to neighbor, RADius of neighbor
	PDCoord  rel_pos;		// vector to neighbor
	for (i = 0; i < nnb; ++i)       // over neighbors
	{	  
#ifdef ATOMVOL
		if (withVol != 0)
		{
			volnb[i] = 0.0;
			fknot[i] = 0;
		}
#endif
		nbs[i]->visitedAs = i;
		rel_pos = nbs[i]->position - pos;       // vector to neighbor
		dist = rel_pos.norm();                  // distance to neighbor
		if (dist == 0)
		{
			std::cerr << "PowerSasa: Invalid distance to neighbour" << std::endl;
			throw PowerSasaException();
		}
		costheta[i] = 1.0;                      // will be overwritten
		nb_RAD = nbs[i]->r;                     // neighbor RADius
		nb_RAD2[i] = nbs[i]->r2;                // neighbor RADius^2
		np[i] = 0;				// initialize number of points
		nt[i] = 0;                              // ?????

		if (dist <= nb_RAD - RAD) // totally covered
		{
			return;          
		}

		if (dist >= RAD + nb_RAD || dist <= RAD - nb_RAD)
		{
			++n_apart;
			np[i] = -1;
			continue;         // neighbour i does not contribute
		}

		costheta[i] = 0.5 * ( dist + (RAD2 - nb_RAD2[i]) / dist ) / RAD;

		if (costheta[i] <= -1.0)  // totally covered (impossible but still...)
		{
			return;          
		}
		if (costheta[i] >= 1.0)   // impossible but still...
		{
			++n_apart;
			np[i] = -1;
			continue;         // neighbour i does not contribute
		}

		sintheta[i] = std::sqrt(1.0 - costheta[i] * costheta[i]);
		e[i] = rel_pos / dist;
		nb_dist[i] = dist;
	}

	if (n_apart == nnb) // no contributing neighbours
	{
		if (withSasa != 0) Sasa[iatom] = 4.0 * pi() * RAD2;
		if (withVol != 0)  Vol[iatom] = (4.0/3.0) * pi() * RAD * RAD2;
		return; 
	}
	

//------ register surface vertices ------------------------

        int partner[2], ptn, ptn0, ptn1;
        const typename POWER_DIAGRAM::PowerDiagram<PDFloat,PDCoord,3>::zeroPoint* zp;
	PDCoord zp_pos;

	int nvx = 0;
	for (unsigned int k = 0; k < atom.myZeroPoints.size(); ++k)
		if(power_diagram->get_zeroPoints()[atom.myZeroPoints[k]].isValid())
	{
		zp = &(power_diagram->get_zeroPoints()[atom.myZeroPoints[k]]);
		zp_pos = zp->getPos();
		ptn = 0;
		for (int kg = 0; kg < 3; ++kg)
		{
			if (zp->generators[kg] != &atom)
			{
				partner[ptn] = zp->generators[kg]->visitedAs;
				++ptn;
			}
		}
		ptn0 = partner[0];
		ptn1 = partner[1];

		if (zp->pos < 0.0 || 1.0 < zp->pos)
		{
			++nt[ptn0];
			++nt[ptn1];
			continue;
		}

		if ((np[ptn0] < 0) || (np[ptn1] < 0))
		{
			std::cerr << "PowerSasa: Invalid surface vertex" << std::endl;
			throw PowerSasaException();
		}
		if ((np[ptn0] >= static_cast<int>(rang.size())) || (np[ptn1] >= static_cast<int>(rang.size())))
		{
			Resize_PNT( (np[ptn0] > np[ptn1]) ? np[ptn0]+1 : np[ptn1]+1 );
			//std::cerr << "PowerSasa: Number of surface vertices for single neighbor exeeds MAX_PNT=" << MAX_PNT << std::endl;
			//throw PowerSasaException();
		}
		if (nvx >= static_cast<int>(vx.size()))
		{
		        Resize_VX(nvx+1);
			//std::cerr << "PowerSasa: number of surface vertices exeeds MAX_VX=" << MAX_VX << std::endl;
			//throw PowerSasaException();
		}

		vx[nvx] = (zp_pos - pos) / RAD;               // faster
		//vx[nvx] = (zp_pos - pos).normalized();      // more accurate?
		p[ptn0][np[ptn0]] = p[ptn1][np[ptn1]] = nvx;  // s-vertex number

		br_c[nvx][0] = ptn0;
		br_c[nvx][1] = ptn1;
		br_p[nvx][0] = np[ptn0];
		br_p[nvx][1] = np[ptn1];
		++nvx;
		++np[ptn0];
		++np[ptn1];

#ifdef ATOMVOL
		typename POWER_DIAGRAM::PowerDiagram<PDFloat,PDCoord,3>::vertex *node1 = zp->from;
		typename POWER_DIAGRAM::PowerDiagram<PDFloat,PDCoord,3>::vertex *node2 = node1->endPoints[zp->branch];
		if (node1->powerValue < 0.0 && node2->powerValue > 0.0)
		{
			if (fknot[ptn0] == 0)
			{
				fknot[ptn0] = 1;
				knot[ptn0] = node1->position;
			}
			else volnb[ptn0] +=
			  std::abs((node1->position - knot[ptn0]).cross(zp_pos - knot[ptn0]).dot(e[ptn0]));
			if (fknot[ptn1] == 0)
			{
				fknot[ptn1] = 1;
				knot[ptn1] = node1->position;
			}
			else volnb[ptn1] +=
			  std::abs((node1->position - knot[ptn1]).cross(zp_pos - knot[ptn1]).dot(e[ptn1]));
		}
		else if (node1->powerValue > 0.0 && node2->powerValue < 0.0)
		{
			if (fknot[ptn0] == 0)
			{
				fknot[ptn0] = 1;
				knot[ptn0] = node2->position;
			}
			else volnb[ptn0] +=
			  std::abs((node2->position - knot[ptn0]).cross(zp_pos - knot[ptn0]).dot(e[ptn0]));
			if (fknot[ptn1] == 0)
			{
				fknot[ptn1] = 1;
				knot[ptn1] = node2->position;
			}
			else volnb[ptn1] +=
			  std::abs((node2->position - knot[ptn1]).cross(zp_pos - knot[ptn1]).dot(e[ptn1]));
		}
		else if (node1->powerValue > 0.0 && node2->powerValue > 0.0)
		{
			PDFloat dpos = node1->powerValue*(1.0 - zp->pos) /
			  (node2->powerValue*zp->pos + node1->powerValue*(1.0 - zp->pos)) - zp->pos;
			if (fknot[ptn0] == 0)
			{
				fknot[ptn0] = 1;
				knot[ptn0] = zp_pos;
			}
			else volnb[ptn0] += 0.5 *
			  std::abs(dpos*(zp_pos - knot[ptn0]).cross(node2->position - node1->position).dot(e[ptn0]));
			if (fknot[ptn1] == 0)
			{
				fknot[ptn1] = 1;
				knot[ptn1] = zp_pos;
			}
			else volnb[ptn1] += 0.5 *
			  std::abs(dpos*(zp_pos - knot[ptn1]).cross(node2->position - node1->position).dot(e[ptn1]));
		}
		else
		{
			std::cerr << "PowerSasa: Impossible zeroPoint" << std::endl;
			throw PowerSasaException();
		}
#endif
	}

	std::vector<typename POWER_DIAGRAM::PowerDiagram<PDFloat,PDCoord,3>::vertex*> const &nodes = atom.myVertices;
	typename POWER_DIAGRAM::PowerDiagram<PDFloat,PDCoord,3>::vertex *node1, *node2;
	s_boost::array<typename POWER_DIAGRAM::PowerDiagram<PDFloat, PDCoord,3>::cellPtr,4> gen1, gen2;
	typename POWER_DIAGRAM::PowerDiagram<PDFloat,PDCoord,3>::cell *at;

	for (j = 0; j < int(nodes.size()); ++j)
	{
		node1 = nodes[j];
		for (kn = 0; kn < 4; ++kn)
		{
			node2 = node1->endPoints[kn];
			if (node2 > node1) continue;
			if (node1->powerValue > 0.0 || node2->powerValue > 0.0) continue;
			gen1 = node1->generators;
			gen2 = node2->generators;
			if ((gen2[0] != &atom) && (gen2[1] != &atom) &&
			  (gen2[2] != &atom) && (gen2[3] != &atom)) continue;
			ptn = 0;
			for (int kg = 0; kg < 4; ++kg)
			{
				at = gen1[kg];
				if (at != &atom &&
				  (at == gen2[0] || at == gen2[1] || at == gen2[2] || at == gen2[3]))
				{
					partner[ptn] = at->visitedAs;
					++ptn;
				}
			}
			ptn0 = partner[0];
			ptn1 = partner[1];
			++nt[ptn0];
			++nt[ptn1];
#ifdef ATOMVOL
			if (fknot[ptn0] == 0)
			{
				fknot[ptn0] = 1;
				knot[ptn0] = node1->position;
			}
			else volnb[ptn0] +=
			  std::abs((node1->position - knot[ptn0]).cross(node2->position - knot[ptn0]).dot(e[ptn0]));
			if (fknot[ptn1] == 0)
			{
				fknot[ptn1] = 1;
				knot[ptn1] = node1->position;
			}
			else volnb[ptn1] +=
			  std::abs((node1->position - knot[ptn1]).cross(node2->position - knot[ptn1]).dot(e[ptn1]));
#endif
		}
	}

// ----- get angles of s-vertices within each circle ---------

	for (i = 0; i < nnb; i++)
	{
		if (np[i] <= 0) continue;
		if (np[i] % 2 != 0)
		{
			std::cerr << "PowerSasa: odd number of crossing between circles of "<<iatom<<" and "<<atom.neighbours[i]-&atoms[0] << std::endl;
			throw PowerSasaException();
		}
		Get_Ang(np[i], p[i], e[i], sintheta[i], costheta[i], ang[i]);
		Get_Next(np[i], ang[i], next[i], p[i], e[i]);
	}

// ---------- get sasa from contours ---------------------------

        int ic1, ic2, ic_0, ic_1, ip2, ip_next, ivx, count;
	PDFloat phi, co, dirdet, ds1, ds2, scone;
        PDCoord *p_ini, *pt, *pt0, vv;
	PDFloat vol2 = 0.0, vol3 = 0.0;

	PDFloat sasa_ia = 0.0;
	for (int iv = 0; iv < nvx; iv++) off[iv] = 0;

	for (int iv = 0; iv < nvx; iv++)
	{
		if (off[iv]) continue;
		p_ini = &vx[iv];
		ic_0 = br_c[iv][0];
		ic_1 = br_c[iv][1];

		// determine which directoin
		dirdet = (e[ic_1].cross(e[ic_0])).dot(*p_ini);
		if (dirdet == 0.0)
		{
		    std::cerr << "PowerSasa: dirdet == 0.0" << std::endl;
		    throw PowerSasaException();
		}
		//assert(dirdet != 0.0);
		if ( dirdet > 0.0 )
		{
			ic1 = ic_0;
//			ip1 = br_p[iv][0];
			ic2 = ic_1;
			ip2 = br_p[iv][1];
		}
		else
		{
			ic1 = ic_1;
//			ip1 = br_p[iv][1];
			ic2 = ic_0;
			ip2 = br_p[iv][0];
		}

		pt = p_ini;
		if (do_sasa) sasa_ia += 2.0 * pi();
		count = 0;

		do                                 // main loop
		{
			++count;
			if (count > MAX_COUNT)
			{
				std::cerr << "PowerSasa: Wrong contour" << std::endl;
				throw PowerSasaException();
			}
			//assert(count <= MAX_PNT);   // else wrong contour
			ip_next = next[ic2][ip2];
			phi = ang[ic2][ip_next] - ang[ic2][ip2];
			if (phi < 0.0) phi += 2.0 * pi();
			
			if (do_sasa)
			{
				co = (e[ic1].dot(e[ic2]) - costheta[ic1]*costheta[ic2])
					/ (sintheta[ic1] * sintheta[ic2]);

				if (co < -1.0) co = -1.0;
				if (co >  1.0) co = 1.0;
				sasa_ia += phi * costheta[ic2] - acos(co);
			}
			
			off[ p[ic2][ip2] ] = 1;
			ic1 = ic2;

			ivx = p[ic1][ip_next];
			pt0 = pt;
			pt = &vx[ivx];
			
			if (br_c[ivx][0] == ic1)
			{
				ic2 = br_c[ivx][1];
				ip2 = br_p[ivx][1];
			}
			else
			{
				ic2 = br_c[ivx][0];
				ip2 = br_p[ivx][0];
			}

			if (withDSasa != 0) {
				ds1 = 0.5 * RAD * phi * 
					( 1.0 + (nb_RAD2[ic1] - RAD2)/(nb_dist[ic1] * nb_dist[ic1]) );
				ds2 = RAD2 / nb_dist[ic1];
				DSasa_parts[iatom][ic1]+= ds1 * e[ic1] - ds2 * (*pt - *pt0).cross(e[ic1]);
			}

			if (withVol != 0 || withDVol != 0) {
				vv = (pos + RAD * (*pt0)).cross(pos + RAD * (*pt));
				scone = sintheta[ic1]*sintheta[ic1]*(phi - sin(phi));
			}
			
			if (withVol != 0)
			{
				vol2 += costheta[ic1]*scone;
#ifdef ATOMVOL
				if (fknot[ic1] == 0)
				{
				  fknot[ic1] = 1;
				  knot[ic1] = pos+RAD*(*pt0);
				}
				else volnb[ic1] +=
				  std::abs((pos+RAD*(*pt0) - knot[ic1]).cross(pos+RAD*(*pt) - knot[ic1]).dot(e[ic1]));
#else
				vol3 -= pos.dot(vv);
#endif
			}
			
			if (withDVol != 0)
			{
				DVol[iatom]-= 0.5 * (vv + (RAD2*scone) * e[ic1]);
			}


		} while (pt != p_ini);
		
		if (do_sasa && sasa_ia > 4.0*pi()) sasa_ia -= 4.0*pi();
	}

// ---------- sasa from single circles ----------------------------

	PDFloat pw_i, pw_j;
	PDCoord  cc;

	for (i = 0; i < nnb; ++i)
	{
		if (np[i] != 0 || nt[i] != 0) continue;

		ok = 1;
		cc = pos + (RAD*costheta[i])*e[i];
		pw_i = -sintheta[i]*sintheta[i]*RAD2;

		for (j = 0; j < nnb; ++j)
		{
			if (j == i) continue;
			pw_j = (nbs[j]->position - cc).squaredNorm() - nb_RAD2[j];
		        if (pw_j <= pw_i)
			{
				ok = 0;
				break;
			}
		}
		if (ok)
		{
			if (do_sasa)
			{  
				sasa_ia += 2.0 * pi() * (1.0 + costheta[i]);
				if (sasa_ia > 4.0*pi()) sasa_ia -= 4.0*pi();
			}
			if (withDSasa != 0)
			{
				DSasa_parts[iatom][i]=( RAD*pi() * (1.0 + (nb_RAD2[i] - RAD2)/(nb_dist[i]*nb_dist[i])) ) * e[i];
			}
			if (withVol != 0 || withDVol != 0) scone = sintheta[i]*sintheta[i] * 2.0*pi();
			if (withVol != 0)   vol2 += costheta[i] * scone;
			if (withDVol != 0)  DVol[iatom] = DVol[iatom] - (0.5*RAD2*scone) * e[i];
		}

	}

	if (withSasa != 0) Sasa[iatom] = RAD2 * sasa_ia;
	if (withDSasa != 0)
	{
		for (i = 0; i < nnb; ++i)
		{
			DSasa[iatom] -= DSasa_parts[iatom][i];
		}
	}

#ifdef ATOMVOL
	if (withVol != 0)
	{
		for (i = 0; i < nnb; ++i)
		{
			if (fknot[i]  != 0) vol3 += RAD * volnb[i] * costheta[i];
		}
	}
#endif
	if (withVol != 0)   Vol[iatom] = RAD*RAD2*sasa_ia/3.0 + RAD*RAD2*vol2/6.0 + vol3/6.0;

	for (i = 0; i < nnb; ++i)       // over neighbors, set zero again (only necessary if an expansion of power diagram is planned)
	{
		nbs[i]->visitedAs = 0;
	}
	return;
}

//=============================================================


template <class PDFloat, class PDCoord> void PowerSasa<PDFloat, PDCoord>::
calc_sasa_all()
{
	for (unsigned int i = 0; i < power_diagram->get_points().size(); ++i)
	{
		calc_sasa_single(i);
	}
}


}
#endif /* SASA_H_ */
