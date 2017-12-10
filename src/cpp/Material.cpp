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

#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

//	Read material data from stream Input
bool CBarMaterial::Read(ifstream& Input, unsigned int mset)
{
    Input >> nset; // Number of property set

    if (nset != mset + 1)
    {
        cerr << "*** Error *** Material sets must be inputted in order !" << endl
             << "    Expected set : " << mset + 1 << endl
             << "    Provided set : " << nset << endl;

        return false;
    }

    Input >> E >> Area; // Young's modulus and section area

    return true;
}

//	Write material data to Stream
void CBarMaterial::Write(COutputter& output, unsigned int mset)
{
    output << setw(5) << mset + 1 << setw(16) << E << setw(16) << Area << endl;
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
    Input >> nset; // Number of property set

    if (nset != mset + 1)
    {
        cerr << "*** Error *** Material sets must be inputted in order !" << endl
             << "    Expected set : " << mset + 1 << endl
             << "    Provided set : " << nset << endl;

        return false;
    }

    Input >> E >> nu; // Young's modulus and Poisson's ratio

    return true;
}


//	Read material data from stream Input
bool CHexMaterial::Read(ifstream& Input, unsigned int mset)
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
void CHexMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset + 1 << setw(16) << E << setw(16) << nu << endl;
}

void CTriangleMaterial::Write(COutputter& output, unsigned int mset)
{
	output << setw(5) << mset+1 << setw(16) << E << setw(16) << nu << endl;
}

//	Write material data to Stream
void CQuadrilateralMaterial::Write(COutputter& output, unsigned int mset)
{
    output << setw(5) << mset + 1 << setw(16) << E << setw(16) << nu << endl;
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

	return true;
}

bool CPlateMaterial::Read(ifstream& Input, unsigned int mset)
{
    Input >> nset; // Number of property set

    if (nset != mset + 1)
    {
        cerr << "*** Error *** Material sets must be inputted in order !" << endl
             << "    Expected set : " << mset + 1 << endl
             << "    Provided set : " << nset << endl;

        return false;
    }

    Input >> E >> h >> nu; // Young's modulus and height and Poisson's ratio

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

void CPlateMaterial::Write(COutputter& output, unsigned int mset)
{
    output << setw(5) << mset + 1 << setw(16) << E << setw(16) << h << setw(16) << nu << endl;
}


bool CShellMaterial::Read(ifstream& Input, unsigned int mset)
{
    Input >> nset; // Number of property set

    if (nset != mset + 1)
    {
        cerr << "*** Error *** Material sets must be inputted in order !" << endl
             << "    Expected set : " << mset + 1 << endl
             << "    Provided set : " << nset << endl;

        return false;
    }

    Input >> E >> h >> nu; // Young's modulus and height and Poisson's ratio

    return true;
}

void CShellMaterial::Write(COutputter& output, unsigned int mset)
{
    output << setw(5) << mset + 1 << setw(16) << E << setw(16) << h << setw(16) << nu << endl;
}

