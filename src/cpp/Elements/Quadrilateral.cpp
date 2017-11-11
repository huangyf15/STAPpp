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
//  which is 77 BTW
unsigned int CQuadrilateral::SizeOfStiffnessMatrix() { return 12 * (12 + 1) / 2; }

void normalize(double vec[3])
{
    double length = std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
    vec[0] /= length;
    vec[1] /= length;
    vec[2] /= length;
}

void AccumulateEtaPsi(const double& eta, const double& psi, const double& weight, const double* xe,
                      const double* ye, double* ke, const double E, const double v)
{
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
        {-Je[1][0] / DetJe, Je[0][0] / DetJe}  // second row
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

    const double d_33 = (1 - v) / 2.0;
    const double k = E / (1 - v * v) * std::abs(DetJe) * weight;

#ifdef b
#error "macro b is predefined"
#else
#define b(ii, jj) (B[8 * (ii - 1) + (jj - 1)])
#endif
    for (unsigned int j = 1; j <= 8; ++j)
    {
        for (unsigned int i = 1; i <= j; ++i)
        {
            // k_ij = B_ki d_kl B_lj
            ke[i - 1 + (j * (j - 1) / 2)] += k * (b(1, i) * (b(1, j) + v * b(2, j))   // merge 1, v
                                                  + b(2, i) * (b(2, j) + v * b(1, j)) // merge v, 1
                                                  + d_33 * b(3, i) * b(3, j));
        }
    }
#undef b
    return;
}

void convert2d23d(const double* k, double* Matrix, const double* i, const double* j)
{
    // to see how these are generated, see ../../memo/4Q2d33d.nb and 4Q2d23d.py
    Matrix[0] = i[0] * (i[0] * k[0] + j[0] * k[1]) + j[0] * (i[0] * k[1] + j[0] * k[2]);
    Matrix[1] = i[1] * (i[1] * k[0] + j[1] * k[1]) + j[1] * (i[1] * k[1] + j[1] * k[2]);
    Matrix[2] = i[1] * (i[0] * k[0] + j[0] * k[1]) + j[1] * (i[0] * k[1] + j[0] * k[2]);
    Matrix[3] = i[2] * (i[2] * k[0] + j[2] * k[1]) + j[2] * (i[2] * k[1] + j[2] * k[2]);
    Matrix[4] = i[2] * (i[1] * k[0] + j[1] * k[1]) + j[2] * (i[1] * k[1] + j[1] * k[2]);
    Matrix[5] = i[2] * (i[0] * k[0] + j[0] * k[1]) + j[2] * (i[0] * k[1] + j[0] * k[2]);
    Matrix[6] = i[0] * (i[0] * k[5] + j[0] * k[8]) + j[0] * (i[0] * k[8] + j[0] * k[9]);
    Matrix[7] = i[0] * (i[2] * k[3] + j[2] * k[4]) + j[0] * (i[2] * k[6] + j[2] * k[7]);
    Matrix[8] = i[0] * (i[1] * k[3] + j[1] * k[4]) + j[0] * (i[1] * k[6] + j[1] * k[7]);
    Matrix[9] = i[0] * (i[0] * k[3] + j[0] * k[4]) + j[0] * (i[0] * k[6] + j[0] * k[7]);
    Matrix[10] = i[1] * (i[1] * k[5] + j[1] * k[8]) + j[1] * (i[1] * k[8] + j[1] * k[9]);
    Matrix[11] = i[1] * (i[0] * k[5] + j[0] * k[8]) + j[1] * (i[0] * k[8] + j[0] * k[9]);
    Matrix[12] = i[1] * (i[2] * k[3] + j[2] * k[4]) + j[1] * (i[2] * k[6] + j[2] * k[7]);
    Matrix[13] = i[1] * (i[1] * k[3] + j[1] * k[4]) + j[1] * (i[1] * k[6] + j[1] * k[7]);
    Matrix[14] = i[1] * (i[0] * k[3] + j[0] * k[4]) + j[1] * (i[0] * k[6] + j[0] * k[7]);
    Matrix[15] = i[2] * (i[2] * k[5] + j[2] * k[8]) + j[2] * (i[2] * k[8] + j[2] * k[9]);
    Matrix[16] = i[2] * (i[1] * k[5] + j[1] * k[8]) + j[2] * (i[1] * k[8] + j[1] * k[9]);
    Matrix[17] = i[2] * (i[0] * k[5] + j[0] * k[8]) + j[2] * (i[0] * k[8] + j[0] * k[9]);
    Matrix[18] = i[2] * (i[2] * k[3] + j[2] * k[4]) + j[2] * (i[2] * k[6] + j[2] * k[7]);
    Matrix[19] = i[2] * (i[1] * k[3] + j[1] * k[4]) + j[2] * (i[1] * k[6] + j[1] * k[7]);
    Matrix[20] = i[2] * (i[0] * k[3] + j[0] * k[4]) + j[2] * (i[0] * k[6] + j[0] * k[7]);
    Matrix[21] = i[0] * (i[0] * k[14] + j[0] * k[19]) + j[0] * (i[0] * k[19] + j[0] * k[20]);
    Matrix[22] = i[0] * (i[2] * k[12] + j[2] * k[13]) + j[0] * (i[2] * k[17] + j[2] * k[18]);
    Matrix[23] = i[0] * (i[1] * k[12] + j[1] * k[13]) + j[0] * (i[1] * k[17] + j[1] * k[18]);
    Matrix[24] = i[0] * (i[0] * k[12] + j[0] * k[13]) + j[0] * (i[0] * k[17] + j[0] * k[18]);
    Matrix[25] = i[0] * (i[2] * k[10] + j[2] * k[11]) + j[0] * (i[2] * k[15] + j[2] * k[16]);
    Matrix[26] = i[0] * (i[1] * k[10] + j[1] * k[11]) + j[0] * (i[1] * k[15] + j[1] * k[16]);
    Matrix[27] = i[0] * (i[0] * k[10] + j[0] * k[11]) + j[0] * (i[0] * k[15] + j[0] * k[16]);
    Matrix[28] = i[1] * (i[1] * k[14] + j[1] * k[19]) + j[1] * (i[1] * k[19] + j[1] * k[20]);
    Matrix[29] = i[1] * (i[0] * k[14] + j[0] * k[19]) + j[1] * (i[0] * k[19] + j[0] * k[20]);
    Matrix[30] = i[1] * (i[2] * k[12] + j[2] * k[13]) + j[1] * (i[2] * k[17] + j[2] * k[18]);
    Matrix[31] = i[1] * (i[1] * k[12] + j[1] * k[13]) + j[1] * (i[1] * k[17] + j[1] * k[18]);
    Matrix[32] = i[1] * (i[0] * k[12] + j[0] * k[13]) + j[1] * (i[0] * k[17] + j[0] * k[18]);
    Matrix[33] = i[1] * (i[2] * k[10] + j[2] * k[11]) + j[1] * (i[2] * k[15] + j[2] * k[16]);
    Matrix[34] = i[1] * (i[1] * k[10] + j[1] * k[11]) + j[1] * (i[1] * k[15] + j[1] * k[16]);
    Matrix[35] = i[1] * (i[0] * k[10] + j[0] * k[11]) + j[1] * (i[0] * k[15] + j[0] * k[16]);
    Matrix[36] = i[2] * (i[2] * k[14] + j[2] * k[19]) + j[2] * (i[2] * k[19] + j[2] * k[20]);
    Matrix[37] = i[2] * (i[1] * k[14] + j[1] * k[19]) + j[2] * (i[1] * k[19] + j[1] * k[20]);
    Matrix[38] = i[2] * (i[0] * k[14] + j[0] * k[19]) + j[2] * (i[0] * k[19] + j[0] * k[20]);
    Matrix[39] = i[2] * (i[2] * k[12] + j[2] * k[13]) + j[2] * (i[2] * k[17] + j[2] * k[18]);
    Matrix[40] = i[2] * (i[1] * k[12] + j[1] * k[13]) + j[2] * (i[1] * k[17] + j[1] * k[18]);
    Matrix[41] = i[2] * (i[0] * k[12] + j[0] * k[13]) + j[2] * (i[0] * k[17] + j[0] * k[18]);
    Matrix[42] = i[2] * (i[2] * k[10] + j[2] * k[11]) + j[2] * (i[2] * k[15] + j[2] * k[16]);
    Matrix[43] = i[2] * (i[1] * k[10] + j[1] * k[11]) + j[2] * (i[1] * k[15] + j[1] * k[16]);
    Matrix[44] = i[2] * (i[0] * k[10] + j[0] * k[11]) + j[2] * (i[0] * k[15] + j[0] * k[16]);
    Matrix[45] = i[0] * (i[0] * k[27] + j[0] * k[34]) + j[0] * (i[0] * k[34] + j[0] * k[35]);
    Matrix[46] = i[0] * (i[2] * k[25] + j[2] * k[26]) + j[0] * (i[2] * k[32] + j[2] * k[33]);
    Matrix[47] = i[0] * (i[1] * k[25] + j[1] * k[26]) + j[0] * (i[1] * k[32] + j[1] * k[33]);
    Matrix[48] = i[0] * (i[0] * k[25] + j[0] * k[26]) + j[0] * (i[0] * k[32] + j[0] * k[33]);
    Matrix[49] = i[0] * (i[2] * k[23] + j[2] * k[24]) + j[0] * (i[2] * k[30] + j[2] * k[31]);
    Matrix[50] = i[0] * (i[1] * k[23] + j[1] * k[24]) + j[0] * (i[1] * k[30] + j[1] * k[31]);
    Matrix[51] = i[0] * (i[0] * k[23] + j[0] * k[24]) + j[0] * (i[0] * k[30] + j[0] * k[31]);
    Matrix[52] = i[0] * (i[2] * k[21] + j[2] * k[22]) + j[0] * (i[2] * k[28] + j[2] * k[29]);
    Matrix[53] = i[0] * (i[1] * k[21] + j[1] * k[22]) + j[0] * (i[1] * k[28] + j[1] * k[29]);
    Matrix[54] = i[0] * (i[0] * k[21] + j[0] * k[22]) + j[0] * (i[0] * k[28] + j[0] * k[29]);
    Matrix[55] = i[1] * (i[1] * k[27] + j[1] * k[34]) + j[1] * (i[1] * k[34] + j[1] * k[35]);
    Matrix[56] = i[1] * (i[0] * k[27] + j[0] * k[34]) + j[1] * (i[0] * k[34] + j[0] * k[35]);
    Matrix[57] = i[1] * (i[2] * k[25] + j[2] * k[26]) + j[1] * (i[2] * k[32] + j[2] * k[33]);
    Matrix[58] = i[1] * (i[1] * k[25] + j[1] * k[26]) + j[1] * (i[1] * k[32] + j[1] * k[33]);
    Matrix[59] = i[1] * (i[0] * k[25] + j[0] * k[26]) + j[1] * (i[0] * k[32] + j[0] * k[33]);
    Matrix[60] = i[1] * (i[2] * k[23] + j[2] * k[24]) + j[1] * (i[2] * k[30] + j[2] * k[31]);
    Matrix[61] = i[1] * (i[1] * k[23] + j[1] * k[24]) + j[1] * (i[1] * k[30] + j[1] * k[31]);
    Matrix[62] = i[1] * (i[0] * k[23] + j[0] * k[24]) + j[1] * (i[0] * k[30] + j[0] * k[31]);
    Matrix[63] = i[1] * (i[2] * k[21] + j[2] * k[22]) + j[1] * (i[2] * k[28] + j[2] * k[29]);
    Matrix[64] = i[1] * (i[1] * k[21] + j[1] * k[22]) + j[1] * (i[1] * k[28] + j[1] * k[29]);
    Matrix[65] = i[1] * (i[0] * k[21] + j[0] * k[22]) + j[1] * (i[0] * k[28] + j[0] * k[29]);
    Matrix[66] = i[2] * (i[2] * k[27] + j[2] * k[34]) + j[2] * (i[2] * k[34] + j[2] * k[35]);
    Matrix[67] = i[2] * (i[1] * k[27] + j[1] * k[34]) + j[2] * (i[1] * k[34] + j[1] * k[35]);
    Matrix[68] = i[2] * (i[0] * k[27] + j[0] * k[34]) + j[2] * (i[0] * k[34] + j[0] * k[35]);
    Matrix[69] = i[2] * (i[2] * k[25] + j[2] * k[26]) + j[2] * (i[2] * k[32] + j[2] * k[33]);
    Matrix[70] = i[2] * (i[1] * k[25] + j[1] * k[26]) + j[2] * (i[1] * k[32] + j[1] * k[33]);
    Matrix[71] = i[2] * (i[0] * k[25] + j[0] * k[26]) + j[2] * (i[0] * k[32] + j[0] * k[33]);
    Matrix[72] = i[2] * (i[2] * k[23] + j[2] * k[24]) + j[2] * (i[2] * k[30] + j[2] * k[31]);
    Matrix[73] = i[2] * (i[1] * k[23] + j[1] * k[24]) + j[2] * (i[1] * k[30] + j[1] * k[31]);
    Matrix[74] = i[2] * (i[0] * k[23] + j[0] * k[24]) + j[2] * (i[0] * k[30] + j[0] * k[31]);
    Matrix[75] = i[2] * (i[2] * k[21] + j[2] * k[22]) + j[2] * (i[2] * k[28] + j[2] * k[29]);
    Matrix[76] = i[2] * (i[1] * k[21] + j[1] * k[22]) + j[2] * (i[1] * k[28] + j[1] * k[29]);
    Matrix[77] = i[2] * (i[0] * k[21] + j[0] * k[22]) + j[2] * (i[0] * k[28] + j[0] * k[29]);
}

//	Calculate element stiffness matrix
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CQuadrilateral::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());

    const CNode& n1 = *nodes[0];
    const CNode& n2 = *nodes[1];
    const CNode& n3 = *nodes[2];
    const CNode& n4 = *nodes[3];

    // =========================== 3d to 2d ============================
    // make p31 p21
    double p31[3] = {n3.XYZ[0] - n1.XYZ[0], n3.XYZ[1] - n1.XYZ[1], n3.XYZ[2] - n1.XYZ[2]};
    double p21[3] = {n2.XYZ[0] - n1.XYZ[0], n2.XYZ[1] - n1.XYZ[1], n2.XYZ[2] - n1.XYZ[2]};
    // double p32[3] = {n3.XYZ[0] - n2.XYZ[0], n3.XYZ[1] - n2.XYZ[1], n3.XYZ[2] - n2.XYZ[2]};

    // n = p31 cross p21 (normalized)
    double n[3] = {p31[1] * p21[2] - p31[2] * p21[1], p31[2] * p21[0] - p31[0] * p21[2],
                   p31[0] * p21[1] - p31[1] * p21[0]};
    normalize(n);

    // i = normalized p21
    // i is manually set parallel to p21 so that y21 = 0
    double i[3] = {p21[0], p21[1], p21[2]};
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

    // =========================== assembly Ke' =========================
    // generate GN4Q for eta, psi

    double ke[36]; // for ke', size = (8, 8), number of elements = 8*9/2 = 36
    clear(ke, 36);
    double pos = 1 / std::sqrt(3.0f);
    double etas[2] = {-pos, pos};
    double psis[2] = {-pos, pos};
    double weights[2][2] = {{1.0, 1.0}, {1.0, 1.0}};

    CQuadrilateralMaterial* material = dynamic_cast<CQuadrilateralMaterial*>(
        ElementMaterial); // Pointer to material of the element
    const double& E = material->E;
    const double& v = material->nu;
    AccumulateEtaPsi(etas[0], psis[0], weights[0][0], xe, ye, ke, E, v);
    AccumulateEtaPsi(etas[0], psis[1], weights[0][1], xe, ye, ke, E, v);
    AccumulateEtaPsi(etas[1], psis[0], weights[1][0], xe, ye, ke, E, v);
    AccumulateEtaPsi(etas[1], psis[1], weights[1][1], xe, ye, ke, E, v);

    // ======================== assembly Ke (2d to 3d) ======================

    convert2d23d(ke, Matrix, i, j);
    return;
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
