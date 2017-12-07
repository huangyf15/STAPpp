/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.02, October 27, 2017                                        */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Material.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

//	Read material data from stream Input
bool CBarMaterial::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set

	if (nset != mset + 1)
	{
		cerr << "*** Error *** Material sets must be inputted in order !" << endl 
			 << "    Expected set : " << mset + 1 << endl
			 << "    Provided set : " << nset << endl;

		return false;
	}

	Input >> E >> Area;	// Young's modulus and section area

	return true;
}

//	Write material data to Stream
void CBarMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset+1 << setw(16) << E << setw(16) << Area << endl;
}

//	Read material data from stream Input
bool CTriangleMaterial::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set

	if (nset != mset + 1)
	{
		cerr << "*** Error *** Material sets must be inputted in order !" << endl 
			 << "    Expected set : " << mset + 1 << endl
			 << "    Provided set : " << nset << endl;

		return false;
	}

	Input >> E >> nu;	// Young's modulus and Poisson's ratio

	return true;
}

bool CQuadrilateralMaterial::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set

	if (nset != mset + 1)
	{
		cerr << "*** Error *** Material sets must be inputted in order !" << endl 
			 << "    Expected set : " << mset + 1 << endl
			 << "    Provided set : " << nset << endl;

		return false;
	}

	Input >> E >> nu;	// Young's modulus and Poisson's ratio

	return true;
}

//	Write material data to Stream
void CTriangleMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset+1 << setw(16) << E << setw(16) << nu << endl;
}

void CQuadrilateralMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset+1 << setw(16) << E << setw(16) << nu << endl;
}

//	Read material data from stream Input
bool CTimoshenkoMaterial::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set

	if (nset != mset + 1)
	{
		cerr << "*** Error *** Material sets must be inputted in order !" << endl
			<< "    Expected set : " << mset + 1 << endl
			<< "    Provided set : " << nset << endl;

		return false;
	}

	Input >> E >> nu >> Area;	// Young's modulus, Poisson's ratio and Section area
	Input >> Iyy >> Izz; 		// Moment of inertia for bending about local axis
	Input >> Thetay1 >> Thetay2 >> Thetay3;
								// Direction cosine of local axis

	return true;
}

//	Write material data to Stream
void CTimoshenkoMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset + 1;

	// Young's modulus, shear modulus, section area,
	// Moment of inertia for bending about local axis,
	// as well as torsional constant and sectional moment
	output << setw(14) << E << setw(14) << nu << setw(14) << Area << setw(14) 
		   << Iyy << setw(14) << Izz;
	
	// Direction cosine of local axis
	output << setw(14) << Thetay1 << setw(14) << Thetay2 << setw(14) << Thetay3 << endl;
		
}