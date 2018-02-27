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
#ifndef POWER_DIAGRAM_H_
#define POWER_DIAGRAM_H_

#define __power_diagram_internal_timing__ 0
#include "array.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include <algorithm>
#include <limits>
#include <ctime>
//#include "basic_vector_calc.h"


namespace POWER_DIAGRAM
{

class MyException : public std::exception
{
public:
	MyException() {	}
};
class IdenticalPointException : public std::exception {};
class VerticesFullException : public std::exception {};




//return the nths natural number as if "without" would be missing (counting starts with zero)
inline int nth(const int n,const int without)
{
	return n+(without<=n);
}
template <class PDCoord, class PDFloat, class Pos_iterator, class Strength_iterator, const int dimension>
void getBoundingBox(PDCoord& lowestCorner,PDCoord& highestCorner,const unsigned int& size,const Pos_iterator pos_begin,const Strength_iterator  strength_begin,const PDFloat additionalCubeSize=pow(2.0,1.0/dimension)-1)
{
	if(size>0)
	{
		lowestCorner=*pos_begin;
		highestCorner=*pos_begin;

		Pos_iterator pos_end=pos_begin;
		Strength_iterator strength_end=strength_begin;

		for(unsigned int i=0;i<size;i++)
		{
			for(int g=dimension-1;g>=0;g--)
			{
				if((*pos_end)[g]-(*strength_end)<lowestCorner[g])lowestCorner[g]=(*pos_end)[g]-(*strength_end);
				if((*pos_end)[g]+(*strength_end)>highestCorner[g])highestCorner[g]=(*pos_end)[g]+(*strength_end);
			}
			++pos_end;
			++strength_end;
		}


		const PDCoord center=0.5*(lowestCorner+highestCorner);
		lowestCorner+=/*PowerDiagram<PDFloat,PDCoord,dimension>::error((highestCorner+lowestCorner).norm())**/(lowestCorner-center)*additionalCubeSize;
		highestCorner+=/*PowerDiagram<PDFloat,PDCoord,dimension>::error((highestCorner+lowestCorner).norm())**/(highestCorner-center)*additionalCubeSize;
	}
}




struct PowerDiagramRuntimeParams
{
	bool radiiGiven;//otherwise the second input array is interpreted as power
	bool fill_myVertices;//otherwise no cell information is constructed, but only a diagram of vertices
	bool fill_neighbours;//otherwise no neighbourhood-information is extracted from myVertices
	bool fill_zeroPoints;//otherwise no information where power diagram is zero is created
	bool with_warnings;//otherwise powerdiagram is in silent mode and will not tell when reducing powers
	bool without_check;
	PowerDiagramRuntimeParams(const bool _radiiGiven,const bool _fill_myVertices,const bool _fill_neighbours,const bool _fill_zeroPoints,const bool _with_warnings,bool _without_check)
	{
		radiiGiven=_radiiGiven;
		fill_myVertices=_fill_myVertices;
		fill_neighbours=_fill_neighbours;
		fill_zeroPoints=_fill_zeroPoints;
		with_warnings=_with_warnings;
		without_check=_without_check;
	}
};
template <class PDFloat, class PDCoord,class Pos_iterator, class Strength_iterator, class BondTo_iterator>
struct PowerDiagramParams
{
	bool create_vertices;//otherwise the object is only constructed, but no calculation is done
	PowerDiagramRuntimeParams runpar;

	const PDCoord lowestCorner;
	const PDCoord highestCorner;

	const unsigned int size;
	Pos_iterator pos_begin;
	Strength_iterator strength_begin;
	BondTo_iterator bondTo_begin;
	
	PowerDiagramParams with_radiiGiven(const bool yon){runpar.radiiGiven=yon; return *this;};
	PowerDiagramParams with_calculate(const bool yon){create_vertices=yon; return *this;};
	PowerDiagramParams with_myVertices(const bool yon){runpar.fill_myVertices=yon; return *this;};
	PowerDiagramParams with_cells(const bool yon){runpar.fill_neighbours=yon;runpar.fill_myVertices=yon; return *this;};
	PowerDiagramParams with_zeroPoints(const bool yon){runpar.fill_zeroPoints=yon;runpar.fill_myVertices=yon; return *this;};
	PowerDiagramParams with_Warnings(const bool yon){runpar.with_warnings=yon;return *this;};
	PowerDiagramParams without_Check(const bool yon){runpar.without_check=yon;return *this;};
	PowerDiagramParams(const unsigned int size_, Pos_iterator& pos_begin_, Strength_iterator& strength_begin_, BondTo_iterator bondTo_begin_, PDCoord& lc, PDCoord& hc,
				 bool radiiGiven_=1, bool create_vertices_=1, bool fill_cellVertices=1,bool fill_cellNeighbours=1,bool fill_zeroPoints=1,bool withWarnings=1,bool withoutChecks=1):
		create_vertices(create_vertices_),runpar(PowerDiagramRuntimeParams(radiiGiven_,fill_cellVertices,fill_cellNeighbours,fill_zeroPoints,withWarnings,withoutChecks)),lowestCorner(lc),highestCorner(hc),
		size(size_),pos_begin(pos_begin_),strength_begin(strength_begin_),bondTo_begin(bondTo_begin_)
	{
	}
};

template <class kind>
struct twoOf
{
	kind a;
	kind b;
};


template <class PDFloat, class PDCoord,const int dimension>
class PowerDiagram
{
public:
	//Forward declarations for cell and vertex
	struct zeroPoint;
	struct vertex;
	struct cell;
	struct EdgeEnds;
	typedef cell* cellPtr;
	typedef cell const* const_cellPtr;
	typedef vertex* vertexPtr;
	typedef vertex const* const_vertexPtr;
	PDCoord center;
	PDFloat maxr2;
private:
	PowerDiagramRuntimeParams params;
	unsigned int _nVertices;
	unsigned int _nUnused;

	std::vector<vertexPtr> unused;
	std::vector< cell > points;
	std::vector< vertex > vertices;//we use the plane space of a vector because of speed. never push_back or it will cause a realloc !!!
	std::vector< zeroPoint> zeros;//this is a vector of all zero points, but maybe some are inactive (inactives do not appear in myZeros!)
	std::vector<twoOf<int> > circlePairs3D;
	std::vector< cell > sideGenerators;


	unsigned int nRevertVertices;
	unsigned int nRevertZeros;
	unsigned int nRevertPoints;
	int cornerOwners[1<<dimension];

	PDFloat powerErr;
	std::vector<vertexPtr> Replaced;//old vertices that are removed
	std::vector<vertexPtr> Invalids;//old vertices that are removed reloaded
	std::vector<cellPtr> Involved; // all cells which are involved
	std::vector< EdgeEnds> planes;//used for connecting new vertices, stores "open ends"
public :
	struct cell
	{
		int visitedAs;
		PDCoord position;
		PDFloat r;
		PDFloat r2;
		typedef zeroPoint* zeroPointPtr;
		cellPtr bondTo;
		std::vector<cellPtr> neighbours;
		std::vector<vertexPtr> myVertices;
		std::vector<int> myZeroPoints;

		inline cell(PDCoord const& pos, PDFloat const& root,PDFloat const& power,const cellPtr bond):visitedAs(0),position(pos),r(root),r2(power),bondTo(bond)
		{
			myVertices.reserve(12);
		}
		inline cell(PDCoord const& pos,const PDFloat& str,const cellPtr bond):visitedAs(0),position(pos),r(str),r2(str*str),bondTo(bond)
		{
			myVertices.reserve(12);
		}
		inline PDFloat power(const PDCoord& coord) const
		{
			return (position-coord).squaredNorm()-r2;
		}
		inline bool isReal(const PowerDiagram<PDFloat,PDCoord,dimension>& This)
		{
			return (this>=&This.getPoints()[0]&&this<&This.getPoints()[This.getPoints().size()]);
		}
	private:

	};
	inline void AddToInvolved(cell& thit)
	{
		thit.visitedAs=Involved.size();
		Involved.push_back(&thit);
	}
	clock_t t1,t2,t3,t4,t5,t6;
	cell const & get_point(const int n) const { return points[n]; }
	std::vector< cell > const & get_points() const { return points; }
        unsigned int get_point_num ( cell const & my_cell )
        {
            return (&my_cell - &points[0]);
        }
	std::vector< vertex> const & get_vertices() const { return vertices; }

	std::vector< zeroPoint> const & get_zeroPoints() const { return zeros; }


	inline static PDFloat error(const PDFloat &f) 
	{
		if(f>std::numeric_limits<PDFloat>::min()/std::numeric_limits<PDFloat>::epsilon())
			return f*(std::numeric_limits<PDFloat>::epsilon());
		else if(f<-std::numeric_limits<PDFloat>::min()/std::numeric_limits<PDFloat>::epsilon())
			return -f*(std::numeric_limits<PDFloat>::epsilon());
		else return std::numeric_limits<PDFloat>::min()/std::numeric_limits<PDFloat>::epsilon();
	}
template <class Pos_iterator, class Strength_iterator, class BondTo_iterator>
static PowerDiagramParams<PDFloat,PDCoord,Pos_iterator,Strength_iterator,BondTo_iterator> create(unsigned int size, Pos_iterator pos_begin, Strength_iterator strength_begin, BondTo_iterator bondTo_begin)
{
	PDCoord highestCorner;
	PDCoord lowestCorner;
	{
		if(size>=1)
			getBoundingBox<PDCoord,PDFloat,Pos_iterator,Strength_iterator,3>(lowestCorner,highestCorner,size,pos_begin,strength_begin);
		else{std::cout<<"create empty PD(not implemented, yet)"<<std::endl;throw MyException();}
		return PowerDiagramParams<PDFloat,PDCoord,Pos_iterator,Strength_iterator,BondTo_iterator>	(size,pos_begin,strength_begin,bondTo_begin,lowestCorner,highestCorner);
	}
}
template <typename Pos_iterator, typename Strength_iterator, typename BondTo_iterator>
	PowerDiagram(PowerDiagramParams<PDFloat,PDCoord,Pos_iterator,Strength_iterator,BondTo_iterator> _params):center(0.5*(_params.highestCorner+_params.lowestCorner)),params(_params.runpar)
	{
		_nUnused=0;
		nRevertPoints=0;
		nRevertZeros=0;
		nRevertVertices=(1<<dimension);
		Replaced.reserve(64);
		planes.resize(64*64);
		points.reserve(_params.size);
		vertices.reserve(_params.size*32+(1<<dimension));
		vertices.resize(vertices.capacity());

		PDCoord lowest=_params.lowestCorner;
		PDCoord highest=_params.highestCorner;
		buildCube(PDCoord(lowest-center),PDCoord(highest-center));
		//build and connect corners
		
		Pos_iterator pos_it=_params.pos_begin;
		Strength_iterator strength_it=_params.strength_begin;
		BondTo_iterator bondTo_it=_params.bondTo_begin;

		if(_params.size>0)
		{
			if(params.radiiGiven)
			{
				points.push_back(cell(*pos_it-center,*strength_it,NULL));
				for(unsigned int i=1;i<_params.size;i++)
					points.push_back(cell(*(++pos_it)-center,*(++strength_it),&points[*(++bondTo_it)]));
			}
			else
			{
				points.push_back(cell(*pos_it-center,sqrt(*strength_it),*strength_it,NULL));
				for(unsigned int i=1;i<_params.size;i++)
					points.push_back(cell(*(++pos_it)-center,sqrt(*(++strength_it)),*(strength_it),&points[*(++bondTo_it)]));
			}
		}

		if(_params.create_vertices)
			buildVertices(points.size());


		if(params.fill_myVertices)
			FillAllMyVertices();
		if(params.fill_neighbours)
			FillAllNeighbours();
		if(params.fill_zeroPoints)
			FillAllZeroPoints();


	}
	const std::vector<cell >& getPoints()const {return points;}

	void revert()
	{
		//if addmore was used this function can revert to the diagram without the added atoms
		for(typename std::vector<vertex > ::iterator it=vertices.begin()+nRevertVertices;it!=vertices.begin()+_nVertices;++it)
			for(typename s_boost::array<cellPtr,dimension+1>::iterator itg=it->generators.begin()+1;itg!=it->generators.end();++itg)
			if((*itg)->isReal(*this))
				(*itg)->myVertices.pop_back();
		_nVertices=nRevertVertices;
		Involved.clear();
		if(points.size()>nRevertPoints)
			points.erase(points.begin()+nRevertPoints,points.end());
		for(int c=0;c<(1<<dimension);c++)
		{
			vertices[c].generators[0]=cornerOwners[c]+&points[0];
			vertices[c].powerValue=vertices[c].generators[0]->power(vertices[c].position);
		}
		for(typename std::vector<vertexPtr > ::iterator iti=Invalids.begin();iti!=Invalids.end();++iti)
		{
			vertexPtr it=*iti;
				it->invalid=0;
				it->rrv=0;
				for(typename s_boost::array <vertexPtr,dimension+1>::const_iterator it2=it->endPoints.begin()+it->isCorner();it2!=it->endPoints.end();++it2)
				for(int g1=it->isCorner();g1<=dimension;g1++)
				for(int g2=(*it2)->isCorner();g2<=dimension;g2++)
				{
					if(it->generators[nth(0,g1)]==(*it2)->generators[nth(0,g2)]&&it->generators[nth(1,g1)]==(*it2)->generators[nth(1,g2)]&&it->generators[nth(2,g1)]==(*it2)->generators[nth(2,g2)])
					{
//						(*it)->endPoints[g1]=*it2;		//is already set 
						(*it2)->endPoints[g2]=&(*it);
					}
				}
				for(typename s_boost::array<cellPtr,dimension+1>::iterator itg=it->generators.begin();itg!=it->generators.end();++itg)
					if((*itg)->isReal(*this))
					{
						(*itg)->myVertices.push_back(&(*it));
						if((*itg)->visitedAs==0)
						{
							(*itg)->visitedAs=-1;
							Involved.push_back(*itg);
						}
					}
		}
		Invalids.clear();
		if(params.fill_neighbours)
			FillAllNeighboursOfInvolved();

		if(params.fill_zeroPoints)
		{
			for(typename std::vector<cellPtr>::const_iterator it=Involved.begin();it!=Involved.end();++it)
				for(typename std::vector<int>::const_reverse_iterator itz=(*it)->myZeroPoints.rbegin();itz!=(*it)->myZeroPoints.rend();++itz)
					if((*itz)>nRevertZeros)
						(*it)->myZeroPoints.pop_back();
					else break;
			zeros.erase(zeros.begin()+nRevertZeros,zeros.end());
		}
		nRevertPoints=0;
		nRevertVertices=0;
		for(int c=0;c<(1<<dimension);c++)
			cornerOwners[c]=0;
	}
	template <class Pos_iterator, class Strength_iterator>
	inline void recalculate(const Pos_iterator pos_it,const Strength_iterator strength_it,const unsigned int size)
	//does standard deletion and calculation of Vertices, neighbour information, ...
	{
		clearAllmyVertices();
		clear_interna();
		if(size>points.size())
		{
			const cellPtr old=&points.front();
			points.reserve(size);
			if(old!=&points.front())
				for(int i=1;i<points.size();i++)
				{
					points[i].bondTo=points[i].bondTo-old+&points[0];
				}
		}
		zeros.clear();
		PDCoord lowest,highest;
		getBoundingBox<PDCoord,PDFloat,PDCoord const*,PDFloat const*,dimension>(lowest,highest,size,&(*pos_it),&(*strength_it));
		buildCube(lowest-center,highest-center);

		if(params.radiiGiven)
			{
				for(unsigned int i=0;i<points.size();i++)
				{
					points[i].position=(*(pos_it+i))-center;
					points[i].r=(*(strength_it+i));
					points[i].r2=points[i].r*points[i].r;
					points[i].visitedAs=0;
				}

				for(unsigned int i=points.size();i<size;i++)
				{
					points.push_back(cell((*(pos_it+i))-center,(*(strength_it+i)),&points.back()));
				}

			}
			else
			{
				for(unsigned int i=0;i<points.size();i++)
				{
					points[i].position=(*(pos_it+i))-center;
					points[i].r2=(*(strength_it+i));
					points[i].r=sqrt(points[i].r2);
					points[i].visitedAs=0;
				}
				for(unsigned int i=points.size();i<size;i++)
				{
					points.push_back(cell((*(pos_it+i))-center,sqrt(*(strength_it+i)),*(strength_it+i),&points.back()));
				}
			}

		buildVertices(size);
		if(params.fill_myVertices||params.fill_neighbours)
			FillAllMyVertices();
		if(params.fill_neighbours)
			FillAllNeighbours();
		if(params.fill_zeroPoints)
			FillAllZeroPoints();
	}


	template <class Pos_iterator, class Strength_iterator>
	inline void addMore(const Pos_iterator pos_it,const Strength_iterator strength_it,const int _newSize)
	//does standard deletion and calculation of Vertices, neighbour information, ...
	{
		const unsigned int gap=_newSize<points.size()?points.size()-_newSize:1;
		const unsigned int newSize=_newSize<points.size()?points.size()+1:_newSize;
		nRevertVertices=_nVertices;
		nRevertZeros=zeros.size();
		nRevertPoints=points.size();
		for(int c=0;c<(1<<dimension);c++)
			cornerOwners[c]=vertices[c].generators[0]-&points[0];
		const const_cellPtr oldstorage=&points[0];
		points.reserve(newSize);
		if(params.radiiGiven)
			{
				for(unsigned int i=0;i<newSize-nRevertPoints;i++)
				{
					points.push_back(cell((*(pos_it+i))-center,(*(strength_it+i)),NULL));
					points.back().bondTo=(&points.back())-gap;
				}
			}
			else
			{
				for(unsigned int i=0;i<newSize-nRevertPoints;i++)
				{
					points.push_back(cell((*(pos_it+i))-center,sqrt(*(strength_it+i)),*(strength_it+i),&points.back()));
					points.back().bondTo=(&points.back())-gap;
				}
			}
		{
			PDCoord lowest;
			PDCoord highest;
			PDCoord rebuild(0,0,0);
			getBoundingBox<PDCoord,PDFloat,PDCoord const*,PDFloat const*,dimension>(lowest,highest,newSize-nRevertPoints,&(*pos_it),&(*strength_it),0.0);


			for(int d=0;d<dimension;d++)
			{
				if(vertices.begin()->position[d]+center[d]-lowest[d]>rebuild[d])rebuild[d]=vertices.begin()->position[d]+center[d]-lowest[d];
				if(vertices[(1<<dimension)-1].position[d]+center[d]-highest[d]<rebuild[d])rebuild[d]=-(vertices[(1<<dimension)-1].position[d]+center[d]-highest[d]);
			}

			if(rebuild.squaredNorm()>0)
			{
				clearAllmyVertices();
				for(typename std::vector<vertex > ::iterator it=vertices.begin();it!=vertices.begin()+_nVertices;++it)
				{
					it->invalid=0;
					it->rrv=0;
				}

			{//createlike
					for(typename std::vector< cell >::iterator it=points.begin();it!=points.begin()+nRevertPoints;++it)
							it->bondTo=it->bondTo-oldstorage+&points[0];
				buildCube(vertices.begin()->position-2*rebuild,vertices[(1<<dimension)-1].position+2*rebuild);



				buildVertices(nRevertPoints);

				if(params.fill_myVertices||params.fill_neighbours)
					FillAllMyVertices(0,(1<<dimension));
				if(params.fill_neighbours)
					FillAllNeighbours();
				if(params.fill_zeroPoints)
					FillAllZeroPoints(nRevertZeros);
			}
				nRevertVertices=_nVertices;

			}
			else
			{
				if(&points[0]!= oldstorage)
				{
					for(typename std::vector<vertex > ::iterator it=vertices.begin();it!=vertices.begin()+_nVertices;++it)
						for(typename s_boost::array<cellPtr,dimension+1>::iterator itg=it->generators.begin();itg!=it->generators.end();++itg)
							if(*itg>=oldstorage&&*itg<=oldstorage+nRevertPoints)
							{
								*itg=*itg-oldstorage+&points[0];
							}

					for(typename std::vector< cell >::iterator it=points.begin();it!=points.begin()+nRevertPoints;++it)
					{
							it->bondTo=it->bondTo-oldstorage+&points[0];
							for(typename std::vector<cellPtr>::iterator itn=it->neighbours.begin();itn!=it->neighbours.end();++itn)
								(*itn)+=&points[0]-oldstorage;
					}
				}
			}
		}


		buildVertices(newSize,nRevertPoints);

		if(params.fill_myVertices||params.fill_neighbours)
			FillAllMyVertices(nRevertPoints,nRevertVertices);
		if(params.fill_neighbours)
			FillAllNeighboursOfInvolved();
		if(params.fill_zeroPoints)
			FillAllZeroPoints(nRevertVertices,nRevertZeros);
	}

	void addMore(const PDCoord& pos,const PDFloat& radius,const int near)
	{
		addMore(&pos,&radius,near);
	}

	static void make_inputfile(std::vector<PDCoord> const& position,std::vector<PDFloat> const& power)
	{
		for (unsigned int i = 0; i < position.size(); ++i)
		{
			std::cout << position[i].x() << " " <<  position[i].y() << " " << position[i].z() << " " <<power[i] << std::endl;
		}
	}

	void clearAllmyVertices()
	{
		for(typename std::vector< cell >::iterator it=points.begin();it!=points.end();++it)
		{
			it->myVertices.clear();
			it->myZeroPoints.clear();
		}
	}
	void buildCube(const PDCoord& lowest,const PDCoord& highest)
	{
		_nVertices=1<<dimension;
		sideGenerators.clear();
		for(int i=0;i<2*dimension;i++)
			sideGenerators.push_back(cell(PDCoord(0,0,0),0,NULL));
		PDCoord lhc=lowest;
		vertices[0].setTo(lowest);
		for(int j=dimension-1;j>=0;j--)
			vertices[0].generators[j+1]=&sideGenerators[j];
		for(int i=0;i<(1<<dimension);i++)
			vertices[i].rrv=0;
		for(int i=0;i<(1<<dimension);i++)
			vertices[i].invalid=0;
		for(int i=1;i<(1<<dimension);i++)
		{
			int j=0;
			while(lhc[j]==highest[j])
			{
				lhc[j]=lowest[j];
				j++;
			}
			lhc[j]=highest[j];
			vertices[i].setTo(lhc);
			for(j=dimension-1;j>=0;j--)
				vertices[i].generators[j+1]=((lhc[j]==lowest[j])?&sideGenerators[j]:&sideGenerators[j+dimension]);
		}
		for(int i=0;i<(1<<dimension);i++)
			for(int d=0;d<dimension;d++)
			{
				const int ii=i;
				const int j=(ii>>d)%2?ii-(1<<d):ii+(1<<d);
				vertices[i].endPoints[d+1]=&vertices[j];
			}
		for(int i=0;i<(1<<dimension);i++)
		{//:TODO: direct sort or faster?
			for(int g=dimension-1;g>0;g--)
				for(int j=g;j>0;j--)
					if(&vertices[i].generators[j+1]-&vertices[i].generators[j]>0)
					{
						std::swap(vertices[i].generators[j],vertices[i].generators[j+1]);
						std::swap(vertices[i].endPoints[j],vertices[i].endPoints[j+1]);
					}
		}
		for(int d=1;d<=dimension;d++)//n
		{
			vertices[0].generators[d]->myVertices.push_back(&vertices[0]);
			vertices[(1<<dimension)-1].generators[d]->myVertices.push_back(&vertices[0]);
		}
	}
	void buildVertices(const unsigned int& nPoints,const int from=0)
	{
		//	try
		{

			if(points.size()>0)
			{
				maxr2=points[0].r2;
				for(int i=from;i< static_cast<int>(nPoints) ;i++)
				{
					if(maxr2<points[i].r2)
						maxr2=points[i].r2;
				}
				powerErr=1000*error(maxr2);
				if(from==0)
				{
					insertFirst();
				}


				for(unsigned int i=(from==0)?1:from;i<nPoints;i++)
				{
					unsigned int done=1;
					while(1)
					{
						try
						{
							doInsertion(prepareInsertion(points[i]));
							break;
						}catch(const PDFloat& errorScale)
						{if(__power_diagram_internal_timing__){t3+=clock();t4+=clock();}

							cellPtr identicalPoint=NULL;
							done++;
							if(done>100)
								throw MyException();
							//is this an Identical point problem?
							{
								PDFloat mindist=(Involved[1]->position-Involved.front()->position).squaredNorm();
								cellPtr closest=Involved[1];
								for(typename std::vector<cellPtr>::const_iterator it=Involved.begin()+2;it!=Involved.end();++it)
									if(((*it)->position-Involved.front()->position).squaredNorm()<mindist)
									{
										mindist=((*it)->position-Involved.front()->position).squaredNorm();
										closest=*it;
									}
								if(error(Involved.front()->r)>sqrt(mindist))
								{
									identicalPoint=closest;
									if(params.with_warnings)
									{
										std::cout<<"numerical similar point to "<<closest-&points[0]+1<<" found. ";
										std::cout<<Involved.front()-&points[0]+1<<" is ignored"<<std::endl;
									}
								}
							}
							//delete new vertices built directly into vertices (when unused has been empty)
							if(_nUnused==0)
							{
								//here comes the deletion ;)
								_nVertices-=Involved.front()->myVertices.size()-unused.size()+_nUnused;
								for(typename std::vector<vertexPtr>::iterator it=Involved.front()->myVertices.begin();it!=Involved.front()->myVertices.end();++it)
									if((*it)-&vertices[0]<(1<<dimension))
									{
										_nVertices++;
									}else (*it)->disconnect();
							}

							//delete new vertices built on unused (not needed because endpoints not constructed,yet)
							for(typename std::vector<vertexPtr>::iterator it=unused.begin()+_nUnused;it!=unused.end();++it)
								(*it)->disconnect();
							_nUnused=unused.size();	

							//reconnect replaced with persisting
							for(typename std::vector<vertexPtr>::const_iterator it=Replaced.begin();it!=Replaced.end();++it)
							for(typename s_boost::array<vertexPtr,dimension+1>::const_iterator it2=(*it)->endPoints.begin()+(*it)->isCorner();it2!=(*it)->endPoints.end();++it2)
			 				if((*it2)->rrv<=0)
							{
								(*it2)->rrv=0;
								for(int g1=(*it)->isCorner();g1<=dimension;g1++)
								for(int g2=(*it2)->isCorner();g2<=dimension;g2++)
								{
									if((*it)->generators[nth(0,g1)]==(*it2)->generators[nth(0,g2)]&&(*it)->generators[nth(1,g1)]==(*it2)->generators[nth(1,g2)]&&(*it)->generators[nth(2,g1)]==(*it2)->generators[nth(2,g2)])	
									{
										(*it)->endPoints[g1]=*it2;
										(*it2)->endPoints[g2]=*it;
									}
								}
							}


							for(typename std::vector<cellPtr>::iterator it=Involved.begin()+1;it!=Involved.end();++it)
								if(!(*it)->myVertices[0]->isConnected())
									(*it)->myVertices[0]=Replaced[0];


							//set everything zero again
							SetInvolvedPersistingVisitedToZero();
							Involved.front()->myVertices.clear();
						 	for(typename std::vector<vertexPtr>::iterator it=Replaced.begin();it!=Replaced.end();++it)
							{
			  					(*it)->rrv=0;
								(*it)->invalid=0;
							}

							Replaced.clear();
							if(identicalPoint!=NULL)
							{	
								Involved.front()->myVertices.push_back(identicalPoint->myVertices.front());
								break;
							}
							else
							{
								const PDFloat oldr2=Involved.front()->r2;
								if(params.with_warnings)
									std::cout<<" Numerical Zero Warning: Power of "<<Involved.front()-&points[0]+1<<" is reduced from "<<Involved.front()->r2;
								Involved.front()->r2-=pow(2.0,done)*(errorScale);
								if(Involved.front()->r2>=0)	Involved.front()->r= sqrt(Involved.front()->r2);
								else 		Involved.front()->r=-sqrt(-Involved.front()->r2);
								if(params.with_warnings)
									std::cout<<" to "<<Involved.front()->r2<<" ( Change was "<<Involved.front()->r2-oldr2<<" )"<<std::endl;
							}




							if(done>100){std::cout<<"Exception : cannot get stable results"<<std::endl; throw MyException();}
						}
					}
				}









if(!params.without_check){
				PDFloat checkconst=0;
//				std::cout<<"checking diagram"<<std::endl;
				for(typename std::vector<vertex > ::iterator it=vertices.begin();it!=vertices.begin()+_nVertices;++it)
					if(it->isConnected())
						for(typename std::vector<cell >::iterator itp=points.begin();itp!=points.end();++itp)
							if(itp->power(it->position)-it->powerValue<-checkconst)
			if((&(*itp)!=it->generators[0])&&(&(*itp)!=it->generators[1])&&(&(*itp)!=it->generators[2])&&(&(*itp)!=it->generators[3]))
					{checkconst=itp->power(it->position)-it->powerValue;
		
						std::cout<<"totaly wrong are "<<itp-points.begin()<<" "<<itp->power(it->position)<<" "<<it->generators[0]->power(it->position)<<" "<<it->generators[1]->power(it->position)<<"          "<<it->generators[0]->position[0]<<" "<<it->generators[0]->position[1]<<" "<<it->generators[0]->position[2]<<std::endl;
std::cout<<it->generators[0]-&points[0]<<" "<<it->generators[1]-&points[0]<<" "<<it->generators[2]-&points[0]<<" "<<it->generators[3]-&points[0]<<std::endl;
std::cout<<itp->power(it->endPoints[0]->position)<<" "<<itp->power(it->endPoints[1]->position)<<std::endl;
std::cout<<it->generators[0]->power(it->endPoints[0]->position)<<" "<<it->generators[0]->power(it->endPoints[1]->position)<<std::endl;
std::cout<<it->generators[1]->power(it->endPoints[0]->position)<<" "<<it->generators[1]->power(it->endPoints[1]->position)<<std::endl;
std::cout<<it->position<<std::endl<<std::endl;
std::cout<<it->endPoints[0]->position<<std::endl<<std::endl;
std::cout<<it->endPoints[1]->position<<std::endl;
exit(1);

}

if(std::abs(checkconst)>0.001)
	std::cout<<"the error of the worst vertex is around "<<checkconst<<std::endl;
//dump_vertices();
}










			}
			else{/*no vertices*/}
		}

		//delete waste that was produced, but during cleanup, never move waste (unvalidating "unused" pointer)
		std::sort(unused.begin(),unused.end());
		for(typename std::vector<vertexPtr>::const_reverse_iterator it=unused.rbegin();it!=unused.rend();++it)
		{
			if((*it)!=&vertices[_nVertices-1])
				vertices[--_nVertices].moveAddressNetworkUpdateOnly(*it);
			else {--_nVertices;}
		}

		unused.clear();
		_nUnused=0;

	}
	inline void findReplacedVertex(vertexPtr& This, PDFloat& value,const cell& insertionPoint)
	{
		if(value<0)
			return;
		PDFloat newValue;
		PDFloat smallVal=std::numeric_limits<PDFloat>::max();//something larger than value
		//look at the neighbours
		for(typename s_boost::array<vertexPtr,dimension+1>::const_iterator it=This->endPoints.begin()+This->isCorner();it!=This->endPoints.end();++it)
		{
			//each powerdiff value defines a plane between insertionPoint and current cell (generator0) approach the direction perpendicular to that plane in direction of insertion point !
			newValue=(*it)->powerdiff3D((*it)->generators[0],&insertionPoint);
		 	if(newValue<value)
			{
				value=newValue;
				This=*it;
				if(value<0)return;
				it=This->endPoints.begin()-(!This->isCorner());
			}
			else if(newValue==value)
				smallVal=newValue;
		}
		if((smallVal!=value))
			return;
		smallVal=std::numeric_limits<PDFloat>::max();
		//found value is not definitely the best one
		//try hard to be sure not beeing in a local minimum (second neighbour)
		for(typename s_boost::array<vertexPtr,dimension+1>::iterator itg=This->endPoints.begin()+This->isCorner();itg!=This->endPoints.end();++itg)
			for(typename s_boost::array<vertexPtr,dimension+1>::iterator itg2=(*itg)->endPoints.begin()+(*itg)->isCorner();itg2!=(*itg)->endPoints.end();++itg2)
				if((*itg2)!=This)
				{
					newValue=(*itg2)->powerdiff3D((*itg2)->generators[0],&insertionPoint);
					if(newValue<value)
					{
						value=newValue;
						This=*itg2;
						return findReplacedVertex(This,value,insertionPoint);
					}
					else if(newValue==value)
						smallVal=newValue;
				}
		if((smallVal!=value))
			return;
		smallVal=std::numeric_limits<PDFloat>::max();
		//second was also close... third neighbour...
		for(typename s_boost::array<vertexPtr,dimension+1>::iterator itg=This->endPoints.begin()+This->isCorner();itg!=This->endPoints.end();++itg)
			for(typename s_boost::array<vertexPtr,dimension+1>::iterator itg2=(*itg)->endPoints.begin()+(*itg)->isCorner();itg2!=(*itg)->endPoints.end();++itg2)
				if((*itg2)!=This)
					for(typename s_boost::array<vertexPtr,dimension+1>::iterator itg3=(*itg2)->endPoints.begin()+(*itg2)->isCorner();itg3!=(*itg2)->endPoints.end();++itg3)
						if((*itg3)!=(*itg)&&(*itg3)!=This)
						{
							newValue=(*itg3)->powerdiff3D((*itg3)->generators[0],&insertionPoint);
							if(newValue<value)
							{
								value=(*itg3)->powerdiff3D((*itg3)->generators[0],&insertionPoint);
								This=*itg3;
								return findReplacedVertex(This,value,insertionPoint);
							}else if(newValue==value)
								smallVal=newValue;
						}
		if(smallVal!=value)
			return;
		if(params.with_warnings)
			std::cout<<"warning : program slowed down because of too small accuracy"<<std::endl;
		//...so the numerical problem wants to be tough? A fat lot we care!
		vertexPtr result=NULL;
		for(typename std::vector<vertex > ::iterator it=vertices.begin();it!=vertices.begin()+nVertices();++it)
		if(it->isConnected())
		{
			if(it->powerdiff3D(it->generators[0],&insertionPoint)<value)
			{
				value=(it)->powerdiff3D((it)->generators[0],&insertionPoint);
				This=&(*it);
			}
		}
		if(result!=NULL)
			This=result;

	}
	inline int finiteReplaced(vertex& This,const_cellPtr const& aCell)
	{
		{

			This.rrv=This.powerdiff3D(aCell,This.generators[0]);
			if(This.rrv>powerErr)return 1;
			else if(This.rrv<-powerErr)return 0;
				else
				{
					This.rrv=0;
					throw &This;
				}
		}
	}
	void dump_vertices(std::ostream& out=std::cout)
	{
		std::cout<<"vertices, generators and neighbours "<<std::endl;

		for(typename std::vector<vertex > ::iterator it=vertices.begin();it!=vertices.begin()+_nVertices;++it)

				if(it->isConnected() && (!it->isCorner()))

				{
					out<<(it->generators[0]-&points[0])<<" "<<std::flush;
					out<<it->generators[1]-&points[0]<<" "<<std::flush;
					out<<it->generators[2]-&points[0]<<" "<<std::flush;
					out<<it->generators[3]-&points[0]<<"    "<<std::flush;
					out<<it->position[0]<<" "<<std::flush;
					out<<it->position[1]<<" "<<std::flush;
					out<<it->position[2]<<" "<<std::flush;
					out<<it->powerValue<<"   "<<std::flush;
					out<<it->generators[0]->power(it->position)<<std::endl;

					out<<" "<<it->endPoints[0]->generators[0]-&points[0]<<" "<<it->endPoints[0]->generators[1]-&points[0]<<" "<<it->endPoints[0]->generators[2]-&points[0]<<" "<<it->endPoints[0]->generators[3]-&points[0]<<std::endl;
					out<<" "<<it->endPoints[1]->generators[0]-&points[0]<<" "<<it->endPoints[1]->generators[1]-&points[0]<<" "<<it->endPoints[1]->generators[2]-&points[0]<<" "<<it->endPoints[1]->generators[3]-&points[0]<<std::endl;
					out<<" "<<it->endPoints[2]->generators[0]-&points[0]<<" "<<it->endPoints[2]->generators[1]-&points[0]<<" "<<it->endPoints[2]->generators[2]-&points[0]<<" "<<it->endPoints[2]->generators[3]-&points[0]<<std::endl;
					out<<" "<<it->endPoints[3]->generators[0]-&points[0]<<" "<<it->endPoints[3]->generators[1]-&points[0]<<" "<<it->endPoints[3]->generators[2]-&points[0]<<" "<<it->endPoints[3]->generators[3]-&points[0]<<std::endl;
		 }else if(it->isCorner()){

//out<<(it-vertices.begin())<<" "<<(it->endPoints[1]-&vertices[0])<<" "<<(it->endPoints[2]-&vertices[0])<<" "<<(it->endPoints[3]-&vertices[0])<<std::endl;
//out<<it->endPoints[3]->position<<std::endl;
					out<<(it->generators[0]-&points[0])<<" "<<(it->generators[1]-&points[0])<<" "<<(it->generators[2]-&points[0])<<" "<<(it->generators[3]-&points[0])<<"    "<<it->position[0]<<" "<<it->position[1]<<" "<<it->position[2]<<" "<<it->powerValue<<"   "<<it->generators[0]->power(it->position)<<std::endl;/*
std::cout<<it->generators[0]->position[0]<<" "<<it->generators[0]->position[1]<<" "<<it->generators[0]->position[2]<<"    "<<it->generators[0]->r2<<std::endl;
std::cout<<it->position[0]<<" "<<it->position[1]<<" "<<it->position[2]<<std::endl;
					out<<" "<<it->endPoints[1]->generators[0]-&points[0]<<" "<<it->endPoints[1]->generators[1]-&points[0]<<" "<<it->endPoints[1]->generators[2]-&points[0]<<" "<<it->endPoints[1]->generators[3]-&points[0]<<std::endl;
					out<<" "<<it->endPoints[2]->generators[0]-&points[0]<<" "<<it->endPoints[2]->generators[1]-&points[0]<<" "<<it->endPoints[2]->generators[2]-&points[0]<<" "<<it->endPoints[2]->generators[3]-&points[0]<<std::endl;
					out<<" "<<it->endPoints[3]->generators[0]-&points[0]<<" "<<it->endPoints[3]->generators[1]-&points[0]<<" "<<it->endPoints[3]->generators[2]-&points[0]<<" "<<it->endPoints[3]->generators[3]-&points[0]<<std::endl;*/
		 }else{std::cout<<"outtake"<<std::endl;}
std::cout<<std::endl;
	}

//	inline const int dimension()const{return dimension;}
	inline int nPoints()const{return points.size();}
	inline unsigned int const& nVertices() const {return _nVertices;}

	bool hasVirtualGenerators(const const_vertexPtr& that)const
	{
		if(!(that->generators[dimension]>=&points[0]&&that->generators[dimension]<&points[points.size()]))
			return 1;
		else	return 0;
	}
	int nVirtualGenerators(const const_vertexPtr& that)const
	{
		if(!(that->generators[dimension]>=&points[0]&&that->generators[dimension]<&points[points.size()]))
			if(!(that->generators[dimension]>=&points[0]&&that->generators[dimension]<&points[points.size()]))
				if(that->isCorner())
					return 3;
				else return 2;
			else	return 1;
		else	return 0;
	}
	cell const *findCellInsideCube(const PDCoord& pos,cell const * hint=NULL)
	{
		if(hint==NULL)hint =&points[points.size()/2];
			for(typename std::vector<cellPtr>::const_iterator it=hint->neighbours.begin();it!=hint->neighbours.end();++it)
			{
				if((*it)->power(pos)<hint->power(pos))return findCellInsideCube(pos,*it);
			}
		return hint;
	}
private:
	vertexPtr getRepresentative(const const_cellPtr This)
	{
		if(This->myVertices.front()->isConnected()&&This->myVertices.front()->hasGenerator(This))
		{
			return This->myVertices.front();
		}
		else 
			for(typename std::vector<vertexPtr>::const_iterator it=This->myVertices.begin()+1;it!=This->myVertices.end();++it)
			{
				if((*it)->isConnected()&&(*it)->hasGenerator(This))
					return *it;
			}
		
		if(!(This==&points[0])) return getRepresentative(This->bondTo);
		else 
		{
			return &vertices[0];
		}
	}
	vertexPtr prepareInsertion(cell & This,vertexPtr hint=NULL)
	{
		//try 
			if(__power_diagram_internal_timing__)t2-=clock();
			//there is a power of new cell that is so low, that only one vertex would be replaced. *hint will be the one
			hint=getRepresentative(This.bondTo);
			PDFloat value=hint->powerdiff3D(hint->generators[0],&This);
			{
				//try
				{
					findReplacedVertex(hint,value,This);
				}/*catch(vertexPtr& trouble){std::cout<<" Numerical Warning: Power of "<<&This-&points[0]<<" reduced from "<<Involved.front()->r2;
						This.r2-=(vertices[0].position-vertices[(1<<dim)].position).norm()*(PowerDiagram<PDFloat,PDCoord>::error());
						if(This.r2>0)	This.r= sqrt(This.r2);
						else 		This.r=-sqrt(-This.r2);
						value=This.bondTo->myVertices.front()->powerdiff(This.bondTo->myVertices.front()->generators[0],&This);
						if(params.with_warning)
							std::cout<<" to "<<Involved.front()->r2<<std::endl;
					}*/

			}

			if(__power_diagram_internal_timing__)t2+=clock();

			try
			{
				FillReplacedPersistingAndInvolved(This,hint);


				if(nRevertPoints>0)
				{//:TODO:
					{
					/*	std::cout<<"replaced "<<Replaced.size()<<std::endl;
						for(typename std::vector<vertex > ::iterator itv=vertices.begin();itv!=vertices.begin()+_nVertices;++itv)
						for(typename std::vector<vertexPtr>::const_iterator it=Replaced.begin();it!=Replaced.end();++it)
							if(*it==&(*itv))
							{
							  std::cout<<(*it)->generators[0]-&points[0]<<" "<<(*it)->generators[1]-&points[0]<<" "<<(*it)->generators[2]-&points[0]<<" "<<(*it)->generators[3]-&points[0]<<"    "<<(*it)->position[0]<<" "<<(*it)->position[1]<<" "<<(*it)->position[2]<<std::endl;
								if(!(*it)->isCorner())std::cout<<" "<<(*it)->endPoints[0]->generators[0]-&points[0]<<" "<<(*it)->endPoints[0]->generators[1]-&points[0]<<" "<<(*it)->endPoints[0]->generators[2]-&points[0]<<" "<<(*it)->endPoints[0]->generators[3]-&points[0]<<"    "<<points.front().power((*it)->endPoints[0]->position)<<" "<<points.back().power((*it)->endPoints[0]->position)<<std::endl;
								std::cout<<" "<<(*it)->endPoints[1]->generators[0]-&points[0]<<" "<<(*it)->endPoints[1]->generators[1]-&points[0]<<" "<<(*it)->endPoints[1]->generators[2]-&points[0]<<" "<<(*it)->endPoints[1]->generators[3]-&points[0]<<"    "<<points.front().power((*it)->endPoints[1]->position)<<" "<<points.back().power((*it)->endPoints[1]->position)<<std::endl;
								std::cout<<" "<<(*it)->endPoints[2]->generators[0]-&points[0]<<" "<<(*it)->endPoints[2]->generators[1]-&points[0]<<" "<<(*it)->endPoints[2]->generators[2]-&points[0]<<" "<<(*it)->endPoints[2]->generators[3]-&points[0]<<"    "<<points.front().power((*it)->endPoints[2]->position)<<" "<<points.back().power((*it)->endPoints[2]->position)<<std::endl;
								std::cout<<" "<<(*it)->endPoints[3]->generators[0]-&points[0]<<" "<<(*it)->endPoints[3]->generators[1]-&points[0]<<" "<<(*it)->endPoints[3]->generators[2]-&points[0]<<" "<<(*it)->endPoints[3]->generators[3]-&points[0]<<"    "<<points.front().power((*it)->endPoints[3]->position)<<" "<<points.back().power((*it)->endPoints[3]->position)<<std::endl;
							}*/
					}
				}
			}
			catch(const_vertexPtr& trouble)
			{
				unsigned int done=1;
				while(1)
				{
					try
					{
						if(done!=1)
						{
							PDFloat value=hint->powerdiff3D(hint->generators[0],&This);
							findReplacedVertex(hint,value,This);
							FillReplacedPersistingAndInvolved(This,hint);
							break;
						}else throw trouble;
					}
					catch(const_vertexPtr& trouble)
					{
						const PDFloat oldr2=Involved.front()->r2;
						if(params.with_warnings)
							std::cout<<"Numerical Warning: Power of "<<Involved.front()-&points[0]+1<<" is reduced from "<<Involved.front()->r2;
						SetInvolvedPersistingVisitedToZero();
						This.myVertices.clear();
						for(typename std::vector<vertexPtr>::const_iterator it=Replaced.begin();it!=Replaced.end();++it)
						{
							(*it)->rrv=0;
							for(int g=(*it)->isCorner();g<=dimension;g++)
								(*it)->endPoints[g]->rrv=0;
						}
						Replaced.clear();
						This.r2-=pow(2.0,done)*(PowerDiagram<PDFloat,PDCoord,dimension>::powerErr);
								if(This.r2>0)	This.r= sqrt(This.r2);
						else 		This.r=-sqrt(-This.r2);
						if(params.with_warnings)
							std::cout<<" to "<<Involved.front()->r2<<" ( Change was "<<Involved.front()->r2-oldr2<<" )"<<std::endl;
						done++;
						if(done>100){std::cout<<"exception : cannot get stable results with atom "<<Involved.front()-&points[0]<<" "<<Involved.front()->position+center<<std::endl;
							throw MyException();}
					}
				}
			}
		return hint;
	}

	void doInsertion(const vertexPtr& hint)
	{
//			if(hint!=NULL)
			{
					if(__power_diagram_internal_timing__){const unsigned int zeit=clock();t3-=zeit;t4-=zeit;}
				CreateFiniteVerticesFromReplaced();
					if(__power_diagram_internal_timing__){const unsigned int zeit=clock();t4+=zeit;t5-=zeit;}
				ConnectNewFinitesAmongThemselves3D();
					if(__power_diagram_internal_timing__)t5+=clock();
				UpdateUnused();
				AssignRepresentativeVerticesToCells(hint);
				SetInvolvedPersistingVisitedToZero();
					if(__power_diagram_internal_timing__)t3+=clock();
			}
	}

	void clear_interna()
	{
		Replaced.clear();
		Involved.clear();
//		planes.clear();//should always be clean
	}

	inline void insertFirst()
	{
		clear_interna();
		for(int i=0;i<(1<<dimension);i++)
			vertices[i].setPowerData(&points.front());
		for(int i=0;i<(1<<dimension);i++)
			vertices[i].generators[0]=&points.front();
		for(int i=0;i<(1<<dimension);i++)
			points.front().myVertices.push_back(&vertices[i]);
	}
	void FillAllMyVertices(const int fromPoint=0,const int fromVertex=1<<dimension)
	{
		{
			if(fromPoint>0)
				Involved.clear();
			for(typename std::vector<cell >::iterator it=points.begin()+fromPoint;it!=points.end();++it)
				it->myVertices.clear();


		for(typename std::vector<vertex >::iterator it=vertices.begin()+fromVertex;it!=vertices.begin()+_nVertices;++it)
		if(!(it->invalid))
		{
			if(!(hasVirtualGenerators(&*it)))
				for(typename s_boost::array<cellPtr,dimension+1>::iterator itg=it->generators.begin();itg!=it->generators.end();++itg)
				{
					(*itg)->myVertices.push_back(&(*it));
					if(fromPoint>0)
						if((*itg)->visitedAs==0)
						{
							(*itg)->visitedAs=-1;
							Involved.push_back(*itg);
						}
				}
				else if(!it->isCorner())
					for(typename s_boost::array<cellPtr,dimension+1>::iterator itg=it->generators.begin();itg!=it->generators.end()&&(!((*itg)-&sideGenerators[0]>=0&&(*itg)-&sideGenerators[0]<2*dimension));++itg)
					{
						(*itg)->myVertices.push_back(&(*it));
						if(fromPoint>0)
							if((*itg)->visitedAs==0)
							{
								(*itg)->visitedAs=-1;
								Involved.push_back(*itg);
							}
					}
					else
					{/*dont give corners to sasa code, it doesnt check for it... so we define myVertices as not holding corners!*/
						std::cout<<"wrong internal order, SASA stopped"<<std::endl;
						throw MyException();
					}

		}
		}
	}

	void FillAllNeighbours()
	{int neighbours=0;
		for(typename std::vector<cell>::iterator it=points.begin();it!=points.end();++it)
		{
			it->neighbours.clear();
			it->visitedAs=-1;
		}
		for(typename std::vector<cell >::iterator it=points.begin();it!=points.end();++it)
			for(typename std::vector<vertexPtr>::const_iterator it2=it->myVertices.begin();it2!=it->myVertices.end();++it2)
				if(!(*it2)->isCorner())
					for(int g=dimension;g>=0;g--)
						if((*it2)->generators[g]->isReal(*this))
						{
							if((*it2)->generators[g]->visitedAs<it-points.begin()&&(*it2)->generators[g]!=&(*it))
							{neighbours++;
								it->neighbours.push_back((*it2)->generators[g]);
								(*it2)->generators[g]->visitedAs=it-points.begin();
							}
						}

		for(typename std::vector<cell >::iterator it=points.begin();it!=points.end();++it)
			it->visitedAs=0;
	}

	void FillAllNeighboursOfInvolved()
	{
		sort(Involved.begin(),Involved.end());
		for(typename std::vector<cellPtr>::iterator it=Involved.begin();it!=Involved.end();++it)
		{
			for(typename std::vector<cellPtr>::iterator itn=(*it)->neighbours.begin();itn!=(*it)->neighbours.end();++itn)
			{
				if((*itn)->visitedAs==-1)
				{
					itn=(*it)->neighbours.erase(itn);
					itn--;
				}
			}
		}
		for(typename std::vector<cellPtr>::iterator it=Involved.begin();it!=Involved.end();++it)
			for(typename std::vector<vertexPtr>::const_iterator it2=(*it)->myVertices.begin();it2!=(*it)->myVertices.end();++it2)
			{
				for(int g=dimension;g>=0;g--)
					if((*it2)->generators[g]->isReal(*this))
					{
						if((*it2)->generators[g]->visitedAs!=0&&(*it2)->generators[g]->visitedAs<=(*it)-&points.front()&&(*it2)->generators[g]!=(*it))
						{
							(*it)->neighbours.push_back((*it2)->generators[g]);
							(*it2)->generators[g]->visitedAs=(*it)-&points.front()+1;
						}
					}
			}

		for(typename std::vector<cellPtr>::iterator it=Involved.begin();it!=Involved.end();++it)
			(*it)->visitedAs=0;
	}
	void FillAllZeroPoints(	unsigned int fromVertex=(1<<dimension),const unsigned int fromZero=0)
	{
		zeros.erase(zeros.begin()+fromZero,zeros.end());
		for(typename std::vector<vertex >::const_iterator it=vertices.begin()+fromVertex;it!=vertices.begin()+this->_nVertices;++it)
			if(!(it->invalid))
				if(it->generators[dimension-1]->isReal(*this))
					for(typename s_boost::array <vertexPtr,dimension+1>::const_iterator it2=it->endPoints.begin()+(hasVirtualGenerators(&*it))*3;it2!=it->endPoints.end();++it2)
					{
						if((*it2)-&(*it)>0)
						{
							if(it->powerValue>0)
							{
								const int branch=it2-it->endPoints.begin();
								const PDFloat& v3=(*it2)->powerValue;
								const PDFloat& v2=it->powerValue;
								const PDFloat& v1=it->generators[branch==0]->power(2*it->position-(*it2)->position);
								const PDFloat quot=2*(v1+v3-2*v2);
								//const PDFloat rootsq=sqr(v1-v3)-4*quot*v2;
								const PDFloat rootsq=(v1-v3)*(v1-v3)-4*quot*v2;
								if(rootsq<=0)
									continue;
								if(quot<powerErr)
									if(v1>=0&&v2>=0&&v3>=0)continue;
								const PDFloat rootquot=sqrt(rootsq)/quot;
								const PDFloat min=(v1-v3)/quot;
								const PDFloat sol1=min+rootquot;
								const PDFloat sol2=min-rootquot;
								if(sol1>0&&sol1<1)
									if((*it2)->powerValue>0)
									{
										zeros.push_back(zeroPoint(it->generators[nth(0,branch)],it->generators[nth(1,branch)],it->generators[nth(2,branch)],sol1,&vertices[it-vertices.begin()],branch));
										zeros.push_back(zeroPoint(it->generators[nth(0,branch)],it->generators[nth(1,branch)],it->generators[nth(2,branch)],sol2,&vertices[it-vertices.begin()],branch));
									}
									else
										zeros.push_back(zeroPoint(it->generators[nth(0,branch)],it->generators[nth(1,branch)],it->generators[nth(2,branch)],sol1,&vertices[it-vertices.begin()],branch));
								else if(sol2>0&&sol2<1)
									zeros.push_back(zeroPoint(it->generators[nth(0,branch)],it->generators[nth(1,branch)],it->generators[nth(2,branch)],sol2,&vertices[it-vertices.begin()],branch));
								else
								{//the covered zeros
									zeros.push_back(zeroPoint(it->generators[nth(0,branch)],it->generators[nth(1,branch)],it->generators[nth(2,branch)],sol1,&vertices[it-vertices.begin()],branch));
									zeros.push_back(zeroPoint(it->generators[nth(0,branch)],it->generators[nth(1,branch)],it->generators[nth(2,branch)],sol2,&vertices[it-vertices.begin()],branch));
								}
							}else if((*it2)->powerValue>0)
							{
								const int branch=it2-it->endPoints.begin();
								const PDFloat& v3=(*it2)->powerValue;
								const PDFloat& v2=it->powerValue;
								const PDFloat& v1=it->generators[branch==0]->power(2*it->position-(*it2)->position);
								const PDFloat quot=2*(v1+v3-2*v2);
								//const PDFloat rootsq=sqr(v1-v3)-4*quot*v2;
								const PDFloat rootsq=(v1-v3)*(v1-v3)-4*quot*v2;
								if(rootsq<=0)
									continue;
								if(quot<powerErr)
									if(v1>=0&&v2>=0&&v3>=0)continue;
								const PDFloat rootquot=sqrt(rootsq)/quot;
								const PDFloat min=(v1-v3)/quot;
								const PDFloat sol1=min+rootquot;
								const PDFloat sol2=min-rootquot;
								if(sol1>0&&sol1<1)
									zeros.push_back(zeroPoint(it->generators[nth(0,branch)],it->generators[nth(1,branch)],it->generators[nth(2,branch)],sol1,&vertices[it-vertices.begin()],branch));
								else
									zeros.push_back(zeroPoint(it->generators[nth(0,branch)],it->generators[nth(1,branch)],it->generators[nth(2,branch)],sol2,&vertices[it-vertices.begin()],branch));
							}
						}
					}
                    /*
					else //vertex on edge of cube dont have connections through spheres
					{

					}
                    */

			for(unsigned int i=fromZero;i<zeros.size();i++)
			{
				zeros[i].generators[0]->myZeroPoints.push_back(i);
				zeros[i].generators[1]->myZeroPoints.push_back(i);
				zeros[i].generators[2]->myZeroPoints.push_back(i);
			}
	}

	inline void tryToBuildVertexOnEdge(const const_vertexPtr& This,const int& here)//,const cellPtr s1, const cellPtr s2,const cellPtr s3,const PDCoord& direction);
	{
		//edge between This (replaced and finite) and that defined by generators s1,s2,s3 will get a vertex (of newest,s1,s2,s3)
		{
			if(_nUnused==0)
			{
				if(nVertices()==vertices.capacity())
				{
					throw  here;
				}				
				vertices[_nVertices].endPointsAndPositionOverwrite(This->endPoints[here],This->getPowerPointOnLine2(This->endPoints[here]));
				vertices[++_nVertices-1].Init(This,here,*this);
			}
			else
			{
				unused[_nUnused-1]->endPointsAndPositionOverwrite(This->endPoints[here],This->getPowerPointOnLine2(This->endPoints[here]));
				unused[--_nUnused]->Init(This,here,*this);
			}
		}

	}


	//  void replace_a_vertex(vertex& old_vertex,const Cell& newGenerator);
	//  void checkvertex(vertexIter& myvertex,Cell& newGenerator,std::vector<vertexIter>& replaced,std::vector <vertexIter>&surroundings,const vertexIter former);

	void FillReplacedPersistingAndInvolved(cell& This,vertexPtr start)
	{
		clear_interna();
		Involved.push_back(&This);
		if(finiteReplaced(*start,Involved.front()))
		{
			if(start->isCorner())
			{
				start->cornerToReplacedAndGo(*this);}
			else
				start->finiteToReplacedAndGo(*this);
		}
	}
	void CreateFiniteVerticesFromReplaced()
	{
		//best procedure for new vertices : knowledge : each new (finite) vertex MUST lie on
		//exactly one old EXISTING edge which is NOT disappearing totally
		//all possible edges are the ones coming out our "replaced" vertices
		//so we only try to create if an endPoint of a replaced vertex is not replaced (visitedAs ==-1)
		for(typename std::vector<vertexPtr>::const_iterator it=Replaced.begin();it!=Replaced.end();++it)
		{
			try
			{
				if(!(*it)->isCorner())
					(*it)->template buildIn<0>(this);
				else
					(*it)->template buildIn<1>(this);
			}
			catch(const int g)
			{
				this->ReserveNewVertices();//at least 2^dim+1 new vertices and at most dim+1 vertices to do, so one realloc is enough, no further try/catch
				for(int g2=g;g2>=(*it)->isCorner();g2--)
				{
					if((*it)->endPoints[g2]->rrv<=0)
					{
                            //this->template tryToBuildVertexOnEdge(*it,g2);
                        this->tryToBuildVertexOnEdge(*it,g2);
					}				
				}
			}
		}
	}

	inline void ConnectNewFinitesAmongThemselves3D()
	{
		//we use an InvolvedSize*InvolvedSize-matrix (sparse) and fill in all new edges
		//these are generated by the new Cell and two older Cells
		//we identify vertices to be connected over an new edge by only the two old cells! (new one is everywhere)
		//the earlier dumps itself into the array at spot [a,b] (with a<b) the later simply picks it up!
		if(Involved.size()*Involved.size()>planes.size())
			planes.resize(Involved.size()*Involved.size());

		for(typename std::vector<vertex >::iterator it=vertices.begin();it!=vertices.begin()+(1<<dimension);++it)

			if(it->rrv>0)
			{
				it->setPowerData(Involved.front());
				it->generators[0]=(Involved.front());
			}


		for(typename std::vector<vertexPtr>::const_iterator it=Involved.front()->myVertices.begin();it!=Involved.front()->myVertices.end();++it)
		//if(!(*it)->isCorner())
		{	
			(*it)->registerForConnection3D(this);
		}

	}

	


	// const int GoAlongFace(const vertexIter& former, const vertexIter& current,const vertexIter& finish,const CellPtr& Generator1,const Ptr& Generator2,const int function(const int&))const;
	void ReserveNewVertices()
	{
		std::vector<int> _replaced;
		std::vector<int> _unused;
		std::vector<int> _currentmyVertices;
		std::vector<int> _first;
		const_vertexPtr const oldMemoryPointer=&vertices.front();
		_replaced.reserve(Replaced.size());
		_unused.reserve(unused.size());
		_currentmyVertices.reserve(Involved.front()->myVertices.size());
		_first.reserve(vertices.capacity());

		for(typename std::vector<vertexPtr>::const_iterator it=Replaced.begin();it!=Replaced.end();++it)
			_replaced.push_back(*it-&vertices.front());
		for(typename std::vector<vertexPtr>::const_iterator it=unused.begin();it!=unused.end();++it)
			_unused.push_back(*it-&vertices.front());
		for(typename std::vector<vertexPtr>::const_iterator it=Involved.front()->myVertices.begin();it!=Involved.front()->myVertices.end();++it)
			_currentmyVertices.push_back(*it-&vertices.front());
		for(typename std::vector<cell >::const_iterator it=points.begin();it!=points.begin()+(Involved.front()-&points.front());++it)
			_first.push_back(it->myVertices.front()-&vertices.front());

		vertices.reserve(2*vertices.capacity()+1);

		for(typename std::vector<vertex >::iterator it=vertices.begin();it!=vertices.begin()+_nVertices;++it)
			it->refreshAfterRealloc(oldMemoryPointer+(&*it-&vertices.front()));
		vertices.resize(vertices.capacity());

		for(unsigned int i=0;i<Replaced.size();i++)
			Replaced[i]=&vertices.front()+_replaced[i];
		for(unsigned int i=0;i<unused.size();i++)
			unused[i]=&vertices.front()+_unused[i];
		for(unsigned int i=0;i<Involved.front()->myVertices.size();i++)
			Involved.front()->myVertices[i]=&vertices.front()+_currentmyVertices[i];
		for(int i=0;i<Involved.front()-&points.front();i++)
			points[i].myVertices.front()=&vertices.front()+_first[i];

	}
	void UpdateUnused()
	{
		unused.resize(_nUnused);
		//mark replaced as unused
		for(typename std::vector<vertexPtr>::iterator it=Replaced.begin();it!=Replaced.end();++it)
		{
			if((*it)->isCorner())
			{
				(*it)->rrv=0;
				//Replaced.erase(it);--it;//the corners are always part of diagram	
				*it=Replaced.back();
				Replaced.pop_back();
				--it;
			}
			else
			{
				if((*it)-&vertices.front()<nRevertVertices)
					Invalids.push_back(*it);
				(*it)->disconnect();//vertex has no connection any more.  we delete it later
			}
		}
		if(nRevertVertices==0)
			unused.insert(unused.end(), Replaced.begin(), Replaced.end());

		else
			for(typename std::vector<vertexPtr>::iterator it=Replaced.begin();it!=Replaced.end();++it)
				if((*it)-&vertices[0]>=nRevertVertices)
					unused.push_back(*it);
				else
					for(typename s_boost::array<cellPtr,dimension+1>::const_iterator itg=(*it)->generators.begin();itg!=(*it)->generators.end();++itg)
						for(unsigned int i=0;i<(*itg)->myVertices.size();i++)
							if((*itg)->myVertices[i]==(*it))
							{
								(*itg)->myVertices.erase((*itg)->myVertices.begin()+i);
								break;
							}
			_nUnused=unused.size();
	}
	void AssignRepresentativeVerticesToCells(const vertexPtr& aDefault) const
	{
		//if there are no new vertices, the new cell is covered
		if(Involved.front()->myVertices.empty())
			Involved.front()->myVertices.push_back(aDefault);
		else// we need one existing vertex close to each cell => we give every cell without representativ a new vertex
			for(typename std::vector<cellPtr>::const_iterator it=Involved.begin()+1;it!=Involved.end();++it)
				if(((*it)>&points[0])&&((*it)<&points[points.size()]))
{
				if(!(*it)->myVertices.front()->isConnected())//if representing vertex has been erased
					(*it)->myVertices.front()=(Involved.front()->myVertices.front());//we assign representative of new also to this one
}else{/*std::cout<<"can never happen"<<std::endl; exit(1);*/}

	}

	void SetInvolvedPersistingVisitedToZero() 
	{
		for(typename std::vector<cellPtr>::const_iterator it=Involved.begin()+1;it!=Involved.end();++it)
			(*it)->visitedAs=0;
		for(typename std::vector<vertexPtr>::const_iterator it=Involved.front()->myVertices.begin();it!=Involved.front()->myVertices.end();++it)
		{
			(*it)->rrv=0;
			if(!(*it)->isCorner())
				(*it)->endPoints[0]->rrv=0;
			else
			{
				(*it)->endPoints[1]->rrv=0;
				(*it)->endPoints[2]->rrv=0;
				(*it)->endPoints[3]->rrv=0;
			}
		}
	}

public:
template <class VectorSubtraction>
    inline /*static*/ PDCoord intersectionOfLineAndPlane3D(const VectorSubtraction& direction,const VectorSubtraction& supportVec,const VectorSubtraction& normal,const PDFloat& planeVal)
	{
		const PDFloat tmp=(normal.dot(direction));
//		const PDFloat sqr=normal.squaredNorm();

        if(tmp>PowerDiagram<PDFloat,PDCoord,dimension>::powerErr*std::numeric_limits<PDFloat>::epsilon())
            return direction*((0.5*(planeVal+normal.squaredNorm())-normal.dot(supportVec))/tmp);
        throw PowerDiagram<PDFloat,PDCoord,dimension>::powerErr;
	}


inline const PDCoord getPowerCenterOf2(const cell *const g0,const cell *const g1)
{
	//0.5*(1+(Ra²-Rb²)/dist2)*(a_pos-b_pos)+a_pos
	return (0.5*(1.+(g0->r2-g1->r2)/(g1->position-g0->position).squaredNorm()))*(g1->position-g0->position)+g0->position;
}

inline PDCoord getPowerPointOnLine(const PDCoord& direction,const PDCoord& supportVector,cell const* const& a, cell const* const& b)
{
//	const PDCoord PlaneNormal=(b->position-a->position)/*/(a->position-b->position).norm()*/;
//	const PDFloat PlaneValue=0.5*(PlaneNormal.squaredNorm()+(a->r2-b->r2)/*(a->position-b->position).norm()*/);
	//PlaneNormal and PlaneValue are a factor of (a->position-b->position).norm() too big but they cancel each other out
	return (PowerDiagram<PDFloat,PDCoord,3>::intersectionOfLineAndPlane3D(direction-supportVector,supportVector-a->position,(b->position-a->position),a->r2-b->r2)+supportVector);
}



struct zeroPoint
{
	PDFloat pos;
	vertexPtr from;
	int branch;
	s_boost::array<cellPtr,dimension> generators;

	zeroPoint(const cellPtr& a,const cellPtr& b,const cellPtr& c,const PDFloat& position,const vertexPtr& origin,const int& way):
		pos(position),from(origin),branch(way)
	{
		generators[0]=a;
//		generators[0]->myZeroPoints.push_back(this);
		generators[1]=b;
//		generators[1]->myZeroPoints.push_back(this);
		generators[2]=c;
//		generators[2]->myZeroPoints.push_back(this);

	}
	PDCoord getPos()const
	{
		return (from->endPoints[branch]->position)*pos-from->position*(pos-1);
	}
	bool isValid()const{return ((!from->invalid)&&(!from->endPoints[branch]->invalid));}
};
struct vertex
{
	PDFloat rrv;//relative replace value (power difference)
	bool invalid;
	s_boost::array<cellPtr,dimension+1> generators;
	PDCoord position;
	PDFloat powerValue;
	s_boost::array <vertexPtr,dimension+1> endPoints;


	friend class PowerDiagram <PDFloat, PDCoord,dimension>;
	inline bool isCorner() const {return endPoints[0]==NULL; }
	inline bool isOnEdge(const PowerDiagram<PDFloat,PDCoord,dimension>& This)
	{
		return (!this->generators[dimension-2]->isReal(This));
	}
	inline bool isOnSurface(const PowerDiagram<PDFloat,PDCoord,dimension>& This)
	{
		return (!this->generators[dimension-1]->isReal(This));
	}
//	inline int hasVirtualGenerators()const {return (generators[dimension]->id<0);}
	inline int hasGenerator(const const_cellPtr& that)const {return (generators[0]==that||generators[1]==that||generators[2]==that||generators[3]==that);}
//	inline int isFinite() const { return generators[0]!=NULL; }
	inline void disconnect(){invalid=1;}
	inline int isConnected()const{return !invalid;}

	//  vertex(const vertex& copy);
	inline vertex():invalid(1)/*,generators(dimension+1,NULL),endPoints(dimension+1,NULL)*/ { }

	inline void Init(const const_vertexPtr& This,const int& keep,const PowerDiagram<PDFloat,PDCoord,dimension>& owner)
	{
			this->setPowerData(owner.Involved.front());

			for(int g=dimension;g>0;g--)
				generators[g]=This->generators[g-(g<=keep)];
			generators[0]=owner.Involved.front();
			owner.Involved.front()->myVertices.push_back(this);
			endPoints[0]->fastWhichis(This)=this;

			if(std::abs(powerValue)<owner.powerErr)
			{
				throw owner.powerErr;
			}
	}

inline PDCoord getPowerPointOnLine2(vertex const* const& persist)const
{
//	const PDCoord PlaneNormal=(b->position-a->position)/*/(a->position-b->position).norm()*/;
//	const PDFloat PlaneValue=0.5*(PlaneNormal.squaredNorm()+(a->r2-b->r2)/*(a->position-b->position).norm()*/);
	//PlaneNormal and PlaneValue are a factor of (a->position-b->position).norm() too big but they cancel each other out
//	const PDFloat lowPower=newOne->power(replaced->position);
//	const PDFloat persistPower=newOne->power(persist->position);
	return ((rrv)/((rrv)-(persist->rrv)))*(persist->position-position)+position;
}

	void operator=(const vertex& that)
	{
		generators=that.generators;
		position=that.position;
		powerValue=that.powerValue;
		endPoints=that.endPoints;
		invalid=that.invalid;
		rrv=that.rrv;
	}
private :
	inline void setPowerData(const const_cellPtr& aCell)
	{
		powerValue=(aCell->position-position).squaredNorm()-aCell->r2;
	}
	inline void setTo(const PDCoord pos)
	{
		position=pos;
		for(int g=dimension;g>=0;g--)
		{
			endPoints[g]=NULL;
			generators[g]=NULL;
		}
	}

	inline	PDFloat powerdiff3D(const_cellPtr const& aCell,const_cellPtr const& bCell)const
	{
		//this has best accuracy when vertex is far away and atoms are close. something similar for a far atom is :
		// -bCall->r2+aCell->r2-(closeCell->position-position).squaredNorm()+2.0*((closeCell->position-position).dot(position-farCell->position)
		return -bCell->r2+aCell->r2-(aCell->position-bCell->position).squaredNorm()+2.0*((aCell->position-bCell->position).dot(position-bCell->position));
	}

	template<class PDCalc>
	inline void endPointsAndPositionOverwrite(const vertexPtr& endPoint,const PDCalc& pos)
	{
		endPoints[0]=endPoint;
		rrv=0;
		invalid=0;
		position=pos;
	}

	void refreshAfterRealloc(const vertex*const& copy)
	{
		for(int g=dimension;g>=0;g--)
			if(this->generators[g]!=NULL&&(!this->generators[g]->myVertices.empty())&&this->generators[g]->myVertices.front()==copy)
				this->generators[g]->myVertices.front()=this;
		for(typename s_boost::array<vertexPtr,dimension+1>::iterator it=endPoints.begin();it!=endPoints.end();++it)
		{
				if((*it)!=NULL)
					*it=this+(*it-copy);
		}
	}

	inline void moveAddressNetworkUpdateOnly(const vertexPtr& whereTo)
	{
		if(this->endPoints[dimension]!=NULL)
			(this->endPoints[dimension]->fastWhichis(this))=whereTo;
		for(typename s_boost::array<vertexPtr,dimension+1>::iterator it=endPoints.begin();it!=endPoints.begin()+dimension;++it)
			((*it)->fastWhichis(this))=whereTo;

		*whereTo=*this;
	}
	inline vertexPtr& fastWhichis (const const_vertexPtr& comp)
	{
		for(typename s_boost::array<vertexPtr,dimension+1>::iterator it=endPoints.begin()+dimension;it!=endPoints.begin();--it)
			if(*it==comp)return *it;
		return endPoints[0];
	}
	inline vertexPtr& persistingWhichis3D (const const_vertexPtr& newOne)
	{
		if(generators[2]==newOne->generators[2])
			if(generators[1]==newOne->generators[1])
				return endPoints[0];
			else
				return endPoints[1];
		else
			if(generators[2]!=newOne->generators[3])
				return endPoints[2];
			else
				return endPoints[3];
	}




	void cornerToReplacedAndGo(PowerDiagram<PDFloat,PDCoord,dimension>& owner)
	{
		owner.Replaced.push_back(this);
		for(typename s_boost::array<cellPtr,dimension+1>::const_iterator it=this->generators.begin();it!=this->generators.end();++it)
			if((*it)->visitedAs==0)
				owner.AddToInvolved(*(*it));
		owner.Involved.front()->myVertices.push_back(this);//although replaced it will be part of the new cell!its a corner!


		for(typename s_boost::array<vertexPtr,dimension+1>::const_iterator it=this->endPoints.begin()+dimension;it!=this->endPoints.begin();--it)
			if((*it)->rrv==0)
			{
				(*it)->replaceCheck(owner);
			}else{}

	}
	void finiteToReplacedAndGo(PowerDiagram<PDFloat,PDCoord,dimension>& owner)
	{
		owner.Replaced.push_back(this);
		for(typename s_boost::array<cellPtr,dimension+1>::const_iterator it=this->generators.begin();it!=this->generators.end();++it)
			if((*it)->visitedAs==0)
				owner.AddToInvolved(*(*it));

		for(typename s_boost::array<vertexPtr,dimension+1>::const_iterator it=this->endPoints.begin();it!=this->endPoints.end();++it)
			if((*it)->rrv==0)
				(*it)->replaceCheck(owner);

	}
	inline void replaceCheck( PowerDiagram<PDFloat,PDCoord,dimension>& owner)
	{
		if(this->isCorner())
			this->cornerReplaceCheck(owner);
		else this->finiteReplaceCheck(owner);
	}

	void finiteReplaceCheck(PowerDiagram<PDFloat,PDCoord,dimension>& owner)
	{
		if(owner.finiteReplaced(*this,owner.Involved.front()))
			this->finiteToReplacedAndGo(owner);
		else
		{
		//	visitedAs=-1;
		}
	}
	void cornerReplaceCheck(PowerDiagram<PDFloat,PDCoord,dimension>& owner)
	{
		if(owner.finiteReplaced(*this,owner.Involved.front()))
			this->cornerToReplacedAndGo(owner);
		else
		{
		//	visitedAs=-1;
		}
	}
	template <const int cornerInfo>
	void buildIn(PowerDiagram<PDFloat,PDCoord,dimension>*const& pd)const
	{
		for(int g=dimension;g>=cornerInfo;g--)
		{
			if(this->endPoints[g]->rrv<=0)
			{
				{
					pd->tryToBuildVertexOnEdge(this,g);
				}
			}
		}
	}


	inline void registerForConnection3D(PowerDiagram<PDFloat,PDCoord,dimension>*const& owner)
	{
		owner->planes[generators[2]->visitedAs*owner->Involved.size()+generators[1]->visitedAs].storeOrConnect(this,endPoints[3]);
		owner->planes[generators[3]->visitedAs*owner->Involved.size()+generators[1]->visitedAs].storeOrConnect(this,endPoints[2]);
		owner->planes[generators[3]->visitedAs*owner->Involved.size()+generators[2]->visitedAs].storeOrConnect(this,endPoints[1]);
	}


};
struct EdgeEnds
{
	vertexPtr a;
	vertexPtr* b;
	inline void storeOrConnect(const vertexPtr& pvertex, vertexPtr& itsEndPointStorage)
	{
		if(this->a==NULL)
		{
			this->a=pvertex;								//we store ourself
			this->b=&itsEndPointStorage; //and where the other should write itself into
		}
		else
		{
			itsEndPointStorage=this->a;//we connect ourself to the other
			*(this->b)=pvertex;								//and the other to us
			this->a=NULL;
		}
	}
	inline void connect(const vertexPtr pvertex, vertexPtr& itsEndPointStorage)
	{
		if(this->a!=NULL)
		{
			itsEndPointStorage=this->a;//we connect ourself to the other
			*(this->b)=pvertex;								//and the other to us
			this->a=NULL;
		}
	}
};

};
//static initializations
//template <class PDFloat, class PDCoord>
//std::vector<cell<PDFloat,PDCoord> > PowerDiagram<PDFloat,PDCoord>::sideGenerators; //outside 



} //namespace POWER_DIAGRAM
#endif /* POWER_DIAGRAM_H_ */
