/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.02, October 27, 2017                                        */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Elements/Quadrilateral.h"

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace std;

//	Constructor
CQuadrilateral::CQuadrilateral()
{
    NEN = 4; // Each element has 4 nodes
    nodes = new CNode*[NEN];

    ND = 12; // 12 DOF in total
    LocationMatrix = new unsigned int[ND];

    ElementMaterial = NULL;
}

//	Desconstructor
CQuadrilateral::~CQuadrilateral() {}

//	Read element data from stream Input
bool CQuadrilateral::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets,
                          CNode* NodeList)
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

    unsigned int MSet;           // Material property set number
    unsigned int N1, N2, N3, N4; // node indexs

    Input >> N1 >> N2 >> N3 >> N4 >> MSet;
    ElementMaterial = dynamic_cast<CQuadrilateralMaterial*>(MaterialSets) + MSet - 1;
    nodes[0] = NodeList + N1 - 1;
    nodes[1] = NodeList + N2 - 1;
    nodes[2] = NodeList + N3 - 1;
    nodes[3] = NodeList + N4 - 1;

    return true;
}

//	Write element data to stream
void CQuadrilateral::Write(COutputter& output, unsigned int Ele)
{
    output << setw(5) << Ele + 1               // element index
           << setw(11) << nodes[0]->NodeNumber // node indexs
           << setw(9) << nodes[1]->NodeNumber  //
           << setw(9) << nodes[2]->NodeNumber  //
           << setw(9) << nodes[3]->NodeNumber  //
           << setw(12) << ElementMaterial->nset << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the
//  element
//	Caution:  Equation number is numbered from 1 !
void CQuadrilateral::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN; N++)
        for (unsigned int D = 0; D < 3; D++)
            LocationMatrix[i++] = nodes[N]->bcode[D];
}

//	Return the size of the element stiffness matrix (stored as an array column by column)
unsigned int CQuadrilateral::SizeOfStiffnessMatrix() { return 12 * (12 + 1) / 2; }

//	Calculate element stiffness matrix
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CQuadrilateral::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());

    const CNode& n1 = *nodes[0];
    const CNode& n2 = *nodes[1];
    const CNode& n3 = *nodes[2];
    const CNode& n4 = *nodes[3];

    // ====================== 3d to 2d ============================
    // make p31 p21
    double p31[3] = {n3.XYZ[0] - n1.XYZ[0], n3.XYZ[1] - n1.XYZ[1], n3.XYZ[2] - n1.XYZ[2]};
    double p21[3] = {n2.XYZ[0] - n1.XYZ[0], n2.XYZ[1] - n1.XYZ[1], n2.XYZ[2] - n1.XYZ[2]};
    // double p32[3] = {n3.XYZ[0] - n2.XYZ[0], n3.XYZ[1] - n2.XYZ[1], n3.XYZ[2] - n2.XYZ[2]};

    // n = p31 cross p21 (normalized)
    double const n[3] = {p31[1] * p21[2] - p31[2] * p21[1], p31[2] * p21[0] - p31[0] * p21[2],
                         p31[0] * p21[1] - p31[1] * p21[0]};
    normalize(n);

    // i = normalized p21
    // i is manually set parallel to p21 so that y21 = 0
    double const i[3] = {p21[0], p21[1], p21[1]};
    normalize(i);
    // j = n cross i
    double const j[3] = {n[1] * i[2] - n[2] * i[1], n[2] * i[0] - n[0] * i[2],
                         n[0] * i[1] - n[1] * i[0]};

    // by here, a conversion matrix is formed,
    // as (x', y') = ((i0, i1, i2), (j0, j1, j2)) . (x, y, z)

    // generate xe, ye
    double xe[4] = {
        i[0] * n1.XYZ[0] + i[1] * n1.XYZ[1] + i[2] * n1.XYZ[2],
        i[0] * n2.XYZ[0] + i[1] * n2.XYZ[1] + i[2] * n2.XYZ[2],
        i[0] * n3.XYZ[0] + i[1] * n3.XYZ[1] + i[2] * n3.XYZ[2],
        i[0] * n4.XYZ[0] + i[1] * n4.XYZ[1] + i[2] * n4.XYZ[2],
    };
    double ye[4] = {
        j[0] * n1.XYZ[0] + j[1] * n1.XYZ[1] + j[2] * n1.XYZ[2],
        j[0] * n2.XYZ[0] + j[1] * n2.XYZ[1] + j[2] * n2.XYZ[2],
        j[0] * n3.XYZ[0] + j[1] * n3.XYZ[1] + j[2] * n3.XYZ[2],
        j[0] * n4.XYZ[0] + j[1] * n4.XYZ[1] + j[2] * n4.XYZ[2],
    };

    // ====================== assembly Ke' =========================
    // generate GN4Q for eta, psi

    double GN4Q[8] = {
        (eta - 1) / 4, (1 - eta) / 4,  (1 + eta) / 4, (-eta - 1) / 4, // first row
        (psi - 1) / 4, (-psi - 1) / 4, (1 + psi) / 4, (1 - psi) / 4   // second row
    };
    double Je[2][2] = {
        {GN4Q[0] * xe[0] + GN4Q[1] * xe[1] + GN4Q[2] * xe[2] + GN4Q[3] * xe[3],
         GN4Q[0] * ye[0] + GN4Q[1] * ye[1] + GN4Q[2] * ye[2] + GN4Q[3] * ye[3]}, // first row
        {GN4Q[4] * xe[0] + GN4Q[5] * xe[1] + GN4Q[6] * xe[2] + GN4Q[7] * xe[3],
         GN4Q[4] * ye[0] + GN4Q[5] * ye[1] + GN4Q[6] * ye[2] + GN4Q[7] * ye[3]} // second row
    };
    double DetJe = Je[0][0] * Je[1][1] - Je[0][1] * Je[1][0];
    double JeI[2][2] = {
        {Je[1][1] / DetJe, -Je[0][1] / DetJe}, // first row
        {Je[1][0] / DetJe, Je[0][0] / DetJe}   // second row
    };

    // generate [B] = Je^-1 * GN
    double BB[8] = {
        JeI[0][0] * GN4Q[0] + JeI[0][1] * GN4Q[4], // Be[1, 1]
        JeI[0][0] * GN4Q[1] + JeI[0][1] * GN4Q[5], // Be[1, 2]
        JeI[0][0] * GN4Q[2] + JeI[0][1] * GN4Q[6], // Be[1, 3]
        JeI[0][0] * GN4Q[3] + JeI[0][1] * GN4Q[7], // first row end
        JeI[1][0] * GN4Q[0] + JeI[1][1] * GN4Q[4], // Be[2, 1]
        JeI[1][0] * GN4Q[1] + JeI[1][1] * GN4Q[5], // Be[2, 2]
        JeI[1][0] * GN4Q[2] + JeI[1][1] * GN4Q[6], // Be[2, 3]
        JeI[1][0] * GN4Q[3] + JeI[1][1] * GN4Q[7]  // second row
    };

    // size of B: 3 by 8
    double B[24] = {
        BB[0], 0,     BB[1], 0,     BB[2], 0,     BB[3], 0,     // first row
        0,     BB[4], 0,     BB[5], 0,     BB[6], 0,     BB[7], // second row
        BB[4], BB[0], BB[5], BB[1], BB[6], BB[2], BB[7], BB[3]  // last row
    };

    CTriangleMaterial* material =
        dynamic_cast<CTriangleMaterial*>(ElementMaterial); // Pointer to material of the element

    const double& E = material->E;
    const double& v = material->nu;
    const double d_33 = (1 - v) / 2.0;

    double ke_[36]; // for ke', size = (8, 8), number of elements = 8*9/2 = 36
    clear(ke_, 36);
#ifdef b
#error "macro b is predefined"
#else
#define b(ii, jj) (B[8 * (ii - 1) + (jj - 1)])
    for (unsigned int i = 1; i <= 8; ++i)
    {
        for (unsigned int j = 1; j <= i; ++j)
        {
            // k_ij = m_ki d_kl m_lj
            ke_[i - 1 + (j * (j - 1) / 2)] = (E/(1*v*v)) * 
            (
                b(1, i) * (b(1, j) + v * b(2, j)) + 
                b(2, i) * (b(2, j) + v * b(1, j)) + 
                d_33 * b(3, i) * b(3, j)
            );
        }
    }
#undef b

    // ======================== assembly Ke (2d to 3d) ======================




    //	Calculate bar length
    // double DX[3]; //	dx = x2-x1, dy = y2-y1, dz = z2-z1
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

    // //	Calculate element stiffness matrix

    // CQuadrilateralMaterial* material = dynamic_cast<CQuadrilateralMaterial*>(
    //     ElementMaterial); // Pointer to material of the element

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

//	Calculate element stress
void CQuadrilateral::ElementStress(double* stress, double* Displacement)
{
    // CQuadrilateralMaterial* material = dynamic_cast<CQuadrilateralMaterial*>(
    //     ElementMaterial); // Pointer to material of the element

    // double DX[3];  //	dx = x2-x1, dy = y2-y1, dz = z2-z1
    // double L2 = 0; //	Square of bar length (L^2)

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
