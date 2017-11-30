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

//!	Material class for Timoshenko beam element
class CTimoshenkoMaterial : public CMaterial
{
public:

	double G;       // Shear modulus

	double Area;    // Sectional area of a bar element

	double Iyy;     // Moment of inertia for bending about local y-axis

	double Iyz;     // Moment of inertia for crossing bending

	double Izz;     // Moment of inertia for bending about local z-axis

	double J;       // Torsional constant

	double Gamma;   // Sectional moment

	double Thetay1; // First direction cosine of local y-axis

	double Thetay2; // Second direction cosine of local y-axis

	double Thetay3; // Third direction cosine of local y-axis

	double Thetaz1; // First direction cosine of local z-axis

	double Thetaz2; // Second direction cosine of local z-axis

	double Thetaz3; // Third direction cosine of local z-axis

public:

	//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

	//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};
