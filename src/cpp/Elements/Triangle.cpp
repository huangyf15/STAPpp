/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.02, October 27, 2017                                        */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Elements/Triangle.h"

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace std;

//  Constructor
CTriangle::CTriangle()
{
    NEN = 3; // Each element has 2 nodes
    nodes = new CNode*[NEN];

    ND = 9;
    LocationMatrix = new unsigned int[ND];

    ElementMaterial = NULL;
}

//  Desconstructor
CTriangle::~CTriangle() {}

//  Read element data from stream Input
bool CTriangle::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int N;

    Input >> N; // element number

    if (N != Ele + 1)
    {
        cerr << "*** Error *** Elements must be inputted in order !" << endl
             << "    Expected element : " << Ele + 1 << endl
             << "    Provided element : " << N << endl;

        return false;
    }

    unsigned int MSet;       // Material property set number
    unsigned int N1, N2, N3; // node number

    Input >> N1 >> N2 >> N3 >> MSet;
    ElementMaterial = dynamic_cast<CTriangleMaterial*>(MaterialSets) + MSet - 1;
    nodes[0] = &NodeList[N1 - 1];
    nodes[1] = &NodeList[N2 - 1];
    nodes[2] = &NodeList[N3 - 1];

    return true;
}

//  Write element data to stream
void CTriangle::Write(COutputter& output, unsigned int Ele)
{
    output << setw(5) << Ele + 1 << setw(11) << nodes[0]->NodeNumber << setw(9)
           << nodes[1]->NodeNumber << setw(12) << nodes[2]->NodeNumber << setw(12)
           << ElementMaterial->nset << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the
//  element
//  Caution:  Equation number is numbered from 1 !
void CTriangle::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN; N++)
        for (unsigned int D = 0; D < 3; D++)
            LocationMatrix[i++] = nodes[N]->bcode[D];
}

//  Return the size of the element stiffness matrix (stored as an array column by column)
//  For 2 node bar element, element stiffness is a 6x6 matrix, whose upper triangular part
//  has 21 elements
unsigned int CTriangle::SizeOfStiffnessMatrix() { return 45; }

inline void normalize(double ptr[3])
{
    double sum = std::sqrt((ptr[0] * ptr[0]) + (ptr[1] * ptr[1]) + (ptr[2] + ptr[2]));
    ptr[0] /= sum;
    ptr[1] /= sum;
    ptr[2] /= sum;
}

//  Calculate element stiffness matrix
//  Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CTriangle::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());

    const CNode& n1 = *nodes[0];
    const CNode& n2 = *nodes[1];
    const CNode& n3 = *nodes[2];

    // make p31 p21
    double p31[3] = {n3.XYZ[0] - n1.XYZ[0], n3.XYZ[1] - n1.XYZ[1], n3.XYZ[2] - n1.XYZ[2]};
    double p21[3] = {n2.XYZ[0] - n1.XYZ[0], n2.XYZ[1] - n1.XYZ[1], n2.XYZ[2] - n1.XYZ[2]};

    // n = p31 cross p21 (normalized)
    double const n[3] = {p31[1] * p21[2] - p31[2] * p21[1], p31[2] * p21[0] - p31[0] * p21[2],
                         p31[0] * p21[1] - p31[1] * p21[0]};
    double Aera = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    n[0] /= Aera;
    n[1] /= Aera;
    n[2] /= Aera;
    Aera /= 2.0;
    // i = normalized p21
    normalize(p21);
    const auto& i = p21;
    // j = n cross i
    double const j[3] = {n[1] * i[2] - n[2] * i[1], n[2] * i[0] - n[0] * i[2],
                         n[0] * i[1] - n[1] * i[0]};
    // by here, a conversion matrix is formed,
    // as (x', y') = ((i0, i1, i2), (j0, j1, j2)) . (x, y, z)

    CTriangleMaterial* material =
        dynamic_cast<CTriangleMaterial*>(ElementMaterial); // Pointer to material of the element

    const double& E = material->E;
    const double& v = material->nu;

    double ke[21] = {
        d * x32 * x32,
        0,

    }

    // generate A as <xyz, ij>

    // //  Calculate bar length
    // double DX[3]; //  dx = x2-x1, dy = y2-y1, dz = z2-z1
    // for (unsigned int i = 0; i < 3; i++)
    //     DX[i] = nodes[1]->XYZ[i] - nodes[0]->XYZ[i];

    // double DX2[6]; //  Quadratic polynomial (dx^2, dy^2, dz^2, dx*dy, dy*dz, dx*dz)
    // DX2[0] = DX[0] * DX[0];
    // DX2[1] = DX[1] * DX[1];
    // DX2[2] = DX[2] * DX[2];
    // DX2[3] = DX[0] * DX[1];
    // DX2[4] = DX[1] * DX[2];
    // DX2[5] = DX[0] * DX[2];

    // double L2 = DX2[0] + DX2[1] + DX2[2];
    // double L = sqrt(L2);

    // //  Calculate element stiffness matrix

    // CTriangleMaterial* material =
    //     dynamic_cast<CTriangleMaterial*>(ElementMaterial); // Pointer to material of the element

    // double k = material->E * material->Area / L / L2;

    // Matrix[0] = k * DX2[0];
    // Matrix[1] = k * DX2[1];
    // Matrix[2] = k * DX2[3];
    // Matrix[3] = k * DX2[2];
    // Matrix[4] = k * DX2[4];
    // Matrix[5] = k * DX2[5];
    // Matrix[6] = k * DX2[0];
    // Matrix[7] = -k * DX2[5];
    // Matrix[8] = -k * DX2[3];
    // Matrix[9] = -k * DX2[0];
    // Matrix[10] = k * DX2[1];
    // Matrix[11] = k * DX2[3];
    // Matrix[12] = -k * DX2[4];
    // Matrix[13] = -k * DX2[1];
    // Matrix[14] = -k * DX2[3];
    // Matrix[15] = k * DX2[2];
    // Matrix[16] = k * DX2[4];
    // Matrix[17] = k * DX2[5];
    // Matrix[18] = -k * DX2[2];
    // Matrix[19] = -k * DX2[4];
    // Matrix[20] = -k * DX2[5];
}

//  Calculate element stress
void CTriangle::ElementStress(double* stress, double* Displacement)
{
    // CTriangleMaterial* material =
    //     dynamic_cast<CTriangleMaterial*>(ElementMaterial); // Pointer to material of the element

    // double DX[3];  //  dx = x2-x1, dy = y2-y1, dz = z2-z1
    // double L2 = 0; //  Square of bar length (L^2)

    // for (unsigned int i = 0; i < 3; i++)
    // {
    //     DX[i] = nodes[1]->XYZ[i] - nodes[0]->XYZ[i];
    //     L2 = L2 + DX[i] * DX[i];
    // }

    // double S[6];
    // for (unsigned int i = 0; i < 3; i++)
    // {
    //     S[i] = -DX[i] * material->E / L2;
    //     S[i + 3] = -S[i];
    // }

    // *stress = 0.0;
    // for (unsigned int i = 0; i < 6; i++)
    // {
    //     if (LocationMatrix[i])
    //         *stress += S[i] * Displacement[LocationMatrix[i] - 1];
    // }
}
