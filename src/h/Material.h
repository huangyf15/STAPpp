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

	double nu;	// poisson's ratio

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

//!	Material class for Beam element
class CBeamMaterial : public CMaterial
{
public:

	double nu; // 泊松比
	double a; // 矩形的宽
	double b; // 矩形的长
	double t1;// 右侧厚度
	double t2;// 上方厚度
	double t3;// 左侧厚度
	double t4;// 下方厚度 
	double n1;// Y'轴的x分量
	double n2;// Y'轴的y分量
	double n3;// Y'轴的y分量

public:
	
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input, unsigned int mset);

//!	Write material data to Stream
	virtual void Write(COutputter& output, unsigned int mset);
};
