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
bool CBeamMaterial::Read(ifstream& Input, unsigned int mset)
{
	Input >> nset;	// Number of property set

	if (nset != mset + 1)
	{
		cerr << "*** Error *** Material sets must be inputted in order !" << endl 
			 << "    Expected set : " << mset + 1 << endl
			 << "    Provided set : " << nset << endl;

		return false;
	}

	Input >> E >> nu >> a >> b >> t1 >> t2 >> t3 >> t4 >> n1 >> n2 >> n3;	// 杨氏模量，泊松比，几何参数和Y’轴指向

	return true;
}

//	Write material data to Stream
void CBeamMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset+1 << setw(16) << E << setw(16) << nu << setw(16) << a << setw(16) << b <<setw(16) << t1 <<setw(16) << t2 <<setw(16) << t3 <<setw(16) << t4 <<endl;
}