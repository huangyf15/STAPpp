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

	Input >> E >> G >> Area;	// Young's modulus, shear modulus and Section area
	Input >> Iyy >> Izz >> J; 
								// Moment of inertia for bending about local axis, torsional constant and Sectional moment
	Input >> Thetay1 >> Thetay2 >> Thetay3 >> Thetaz1 >> Thetaz2 >> Thetaz3;
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
	output << setw(16) << E << setw(16) << G << setw(16) << Area << setw(16) 
		   << Iyy << setw(16) << Izz << setw(16) << J << endl;
	
	// Direction cosine of local axis
	output << setw(5) << Thetay1 << setw(16) << Thetay2 << setw(16) << Thetay3 
		<< setw(16) << Thetaz1 << setw(16) << Thetaz2 << setw(16) << Thetaz3 << endl;
		
}