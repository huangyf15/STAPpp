/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.02, October 27, 2017                                        */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "Outputter.h"

#include <stddef.h>
#include <iostream>
#include <fstream>

using namespace std;

//!	Material base class which only define one data member
/*!	All type of material classes should be derived from this base class */
class CMaterial
{
public:

	unsigned int nset;	//!< Number of set
	
	double E;  //!< Young's modulus

public:

//! Virtual deconstructor
    virtual ~CMaterial() {};

//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset) = 0;

//!	Write material data to Stream
    virtual void Write(COutputter& output, unsigned int mset) = 0;

};

//!	Material class for bar element
class CBarMaterial : public CMaterial
{
public:

	double Area;	//!< Sectional area of a bar element

public:
	
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};

//!	Material class for bar element
class CTriangleMaterial : public CMaterial
{
public:

	double nu; // Poisson's ratio

public:
	
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};

//!	Material class for Quadrilateral element
class CQuadrilateralMaterial : public CMaterial
{
public:

	double nu; // Poisson's ratio

public:
	
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};

class CHexMaterial : public CMaterial
{
public:

	double nu;	//!< Sectional area of a bar element

public:

	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

	//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};

//!	Material class for Beam element
class CBeamMaterial : public CMaterial
{
public:

	double nu; // Poisson ratio
	double a; // wide of rectangle
	double b; // height of rectangle
	double t1;// right thickness
	double t2;// above thickness
	double t3;// left thickness
	double t4;// below thickness 
	double n1;// x component of y' axis
	double n2;// y component of y' axis
	double n3;// z component of y' axis

	public:
	
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};

class CPlateMaterial : public CMaterial
{
public:

	double h; //thickness or height

	double nu; // Poisson's ratio

public:
	
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};


//!	Material class for Timoshenko beam element
class CTimoshenkoMaterial : public CMaterial
{
public:

	double nu;      // Poisson's ratio

	double Area;    // Sectional area of a beam element

	double Iyy;     // Moment of inertia for bending about local y-axis

	double Izz;     // Moment of inertia for bending about local z-axis

	double Thetay1; // 1st direction cosine of local y-axis

	double Thetay2; // 2nd direction cosine of local y-axis

	double Thetay3; // 3rd direction cosine of local y-axis

public:

	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

	//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};
