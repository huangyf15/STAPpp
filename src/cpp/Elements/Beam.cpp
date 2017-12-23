/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.02, October 27, 2017                                        */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Elements/Beam.h"

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace std;

//	Constructor
CBeam::CBeam()
{
    NEN = 2; // Each element has 2 nodes
    nodes = new CNode*[NEN];

    ND = 12;
    LocationMatrix = new unsigned int[ND];

    ElementMaterial = nullptr;
}

//	Desconstructor
CBeam::~CBeam() {}

//	Read element data from stream Input
bool CBeam::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
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

    unsigned int MSet;   // Material property set number
    unsigned int N1, N2; // Left node number and right node number
    Input >> N1 >> N2 >> MSet;
    ElementMaterial = static_cast<CBeamMaterial*>(MaterialSets) + MSet - 1;
    nodes[0] = &NodeList[N1 - 1];
    nodes[1] = &NodeList[N2 - 1];

    return true;
}

//	Write element data to stream
void CBeam::Write(COutputter& output, unsigned int Ele)
{
    output << setw(5) << Ele + 1 << setw(11) << nodes[0]->NodeNumber << setw(9)
           << nodes[1]->NodeNumber << setw(12) << ElementMaterial->nset << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the
//  element
//	Caution:  Equation number is numbered from 1 !
void CBeam::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN; N++)
    {
        for (unsigned int D = 0; D < 6; D++)
        {
            LocationMatrix[i++] = nodes[N]->bcode[D];
        }
    }
}

//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For 2 node beam element, element stiffness is a 12x12 matrix, whose upper triangular part
//	has 21 elements
unsigned int CBeam::SizeOfStiffnessMatrix() { return 78; }

//	Calculate element stiffness matrix
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CBeam::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());

    //	Calculate beam length
    double DX[3]; //	dx = x2-x1, dy = y2-y1, dz = z2-z1
    for (unsigned int i = 0; i < 3; i++)
        DX[i] = nodes[1]->XYZ[i] - nodes[0]->XYZ[i];

    double L = sqrt(DX[0] * DX[0] + DX[1] * DX[1] + DX[2] * DX[2]);

    //计算局部坐标系的矩阵需要的元素
    const CBeamMaterial& material =
        static_cast<CBeamMaterial&>(*ElementMaterial); // Pointer to material of the element
    double Iz = (material.a * material.a * material.a) * material.b / 12 -
                pow(material.a - material.t1 - material.t3, 3.0) *
                    (material.b - material.t2 - material.t4) / 12;
    double Iy = (material.b * material.b * material.b) * material.a / 12 -
                pow(material.b - material.t2 - material.t4, 3.0) *
                    (material.a - material.t1 - material.t3) / 12;
    double Ip = Iz + Iy;
    double k1 = material.E * material.a * material.b / L;
    double k2 = 12 * material.E * Iz / (L * L * L);
    double k3 = 12 * material.E * Iy / (L * L * L);
    double k4 = material.E * Ip / ((2 + 2 * material.nu) * L);
    double k5 = 4 * material.E * Iy / L;
    double k6 = 4 * material.E * Iz / L;
    double k7 = 6 * material.E * Iy / (L * L);
    double k8 = 6 * material.E * Iz / (L * L);

    double n[3][3];

    for (unsigned int i = 0; i < 3; i++)
        n[0][i] = DX[i] / L; //orientation of x'

    n[1][0] = material.n1;
    n[1][1] = material.n2;
    n[1][2] = material.n3; // orientation of y'
    n[2][0] = n[0][1] * n[1][2] - n[0][2] * n[1][1];
    n[2][1] = n[0][2] * n[1][0] - n[0][0] * n[1][2];
    n[2][2] = n[0][0] * n[1][1] - n[0][1] * n[1][0]; //orientation of z'

    double N1 = n[0][0] * n[0][0];
    double N2 = n[0][1] * n[0][1];
    double N3 = n[0][2] * n[0][2];
    double N4 = n[1][0] * n[1][0];
    double N5 = n[1][1] * n[1][1];
    double N6 = n[1][2] * n[1][2];
    double N7 = n[2][0] * n[2][0];
    double N8 = n[2][1] * n[2][1];
    double N9 = n[2][2] * n[2][2];
    double N10 = n[0][0] * n[0][1];
    double N11 = n[1][0] * n[1][1];
    double N12 = n[2][0] * n[2][1];
    double N13 = n[0][1] * n[0][2];
    double N14 = n[1][1] * n[1][2];
    double N15 = n[2][1] * n[2][2];
    double N16 = n[0][0] * n[0][2];
    double N17 = n[1][0] * n[1][2];
    double N18 = n[2][0] * n[2][2];
    double N19 = n[1][2] * n[2][0];
    double N20 = n[1][0] * n[2][2];
    double N21 = n[1][1] * n[2][0];
    double N22 = n[1][0] * n[2][1];
    double N23 = n[1][2] * n[2][1];
    double N24 = n[1][1] * n[2][2];
    double N25 = n[1][0] * n[2][0];
    double N26 = n[1][1] * n[2][1];
    double N27 = n[1][2] * n[2][2];

    Matrix[0] = k1 * N1 + k2 * N4 + k3 * N7;
    Matrix[1] = k1 * N2 + k2 * N5 + k3 * N8;
    Matrix[2] = k1 * N10 + k2 * N11 + k3 * N12;
    Matrix[3] = k1 * N3 + k2 * N6 + k3 * N9;
    Matrix[4] = k1 * N13 + k2 * N14 + k3 * N15;
    Matrix[5] = k1 * N16 + k2 * N17 + k3 * N18;
    Matrix[6] = k4 * N1 + k5 * N4 + k6 * N7;
    Matrix[7] = k8 * N19 - k7 * N20;
    Matrix[8] = k8 * N21 - k7 * N22;
    Matrix[9] = k8 * N25 - k7 * N25;
    Matrix[10] = k4 * N2 + k5 * N5 + k6 * N8;
    Matrix[11] = k4 * N10 + k5 * N11 + k6 * N12;
    Matrix[12] = k8 * N23 - k7 * N24;
    Matrix[13] = k8 * N26 - k7 * N26;
    Matrix[14] = k8 * N22 - k7 * N21;
    Matrix[15] = k4 * N3 + k5 * N6 + k6 * N9;
    Matrix[16] = k4 * N13 + k5 * N14 + k6 * N15;
    Matrix[17] = k4 * N16 + k5 * N17 + k6 * N18;
    Matrix[18] = k8 * N27 - k7 * N27;
    Matrix[19] = k8 * N24 - k7 * N23;
    Matrix[20] = k8 * N20 - k7 * N19;
    Matrix[21] = k1 * N1 + k2 * N4 + k3 * N7;
    Matrix[22] = k7 * N19 - k8 * N20;
    Matrix[23] = k7 * N21 - k8 * N22;
    Matrix[24] = k7 * N25 - k8 * N25;
    Matrix[25] = -k1 * N16 - k2 * N17 - k3 * N18;
    Matrix[26] = -k1 * N10 - k2 * N11 - k3 * N12;
    Matrix[27] = -k1 * N1 - k2 * N4 - k3 * N7;
    Matrix[28] = k1 * N2 + k2 * N5 + k3 * N8;
    Matrix[29] = k1 * N10 + k2 * N11 + k3 * N12;
    Matrix[30] = k7 * N23 - k8 * N24;
    Matrix[31] = k7 * N26 - k8 * N26;
    Matrix[32] = k7 * N22 - k8 * N21;
    Matrix[33] = -k1 * N13 - k2 * N14 - k3 * N15;
    Matrix[34] = -k1 * N2 - k2 * N5 - k3 * N8;
    Matrix[35] = -k1 * N10 - k2 * N11 - k3 * N12;
    Matrix[36] = k1 * N3 + k2 * N6 + k3 * N9;
    Matrix[37] = k1 * N13 + k2 * N14 + k3 * N15;
    Matrix[38] = k1 * N16 + k2 * N17 + k3 * N18;
    Matrix[39] = k7 * N27 - k8 * N27;
    Matrix[40] = k7 * N24 - k8 * N23;
    Matrix[41] = k7 * N20 - k8 * N19;
    Matrix[42] = -k1 * N3 - k2 * N6 - k3 * N9;
    Matrix[43] = -k1 * N13 - k2 * N14 - k3 * N15;
    Matrix[44] = -k1 * N16 - k2 * N17 - k3 * N18;
    Matrix[45] = k4 * N1 + k5 * N4 + k6 * N7;
    Matrix[46] = k7 * N20 - k8 * N19;
    Matrix[47] = k7 * N22 - k8 * N21;
    Matrix[48] = k7 * N25 - k8 * N25;
    Matrix[49] = (k5 * N17) / 2 - k4 * N16 + (k6 * N18) / 2;
    Matrix[50] = (k5 * N11) / 2 - k4 * N10 + (k6 * N12) / 2;
    Matrix[51] = -k4 * N1 + (k5 * N4) / 2 + (k6 * N7) / 2;
    Matrix[52] = k8 * N19 - k7 * N20;
    Matrix[53] = k8 * N21 - k7 * N22;
    Matrix[54] = k8 * N25 - k7 * N25;
    Matrix[55] = k4 * N2 + k5 * N5 + k6 * N8;
    Matrix[56] = k4 * N10 + k5 * N11 + k6 * N12;
    Matrix[57] = k7 * N24 - k8 * N23;
    Matrix[58] = k7 * N26 - k8 * N26;
    Matrix[59] = k7 * N21 - k8 * N22;
    Matrix[60] = (k5 * N14) / 2 - k4 * N13 + (k6 * N15) / 2;
    Matrix[61] = -k4 * N2 + (k5 * N5) / 2 + (k6 * N8) / 2;
    Matrix[62] = (k5 * N11) / 2 - k4 * N10 + (k6 * N12) / 2;
    Matrix[63] = k8 * N23 - k7 * N24;
    Matrix[64] = k8 * N26 - k7 * N26;
    Matrix[65] = k8 * N22 - k7 * N21;
    Matrix[66] = k4 * N3 + k5 * N6 + k6 * N9;
    Matrix[67] = k4 * N13 + k5 * N14 + k6 * N15;
    Matrix[68] = k4 * N16 + k5 * N17 + k6 * N18;
    Matrix[69] = k7 * N27 - k8 * N27;
    Matrix[70] = k7 * N23 - k8 * N24;
    Matrix[71] = k7 * N19 - k8 * N20;
    Matrix[72] = -k4 * N3 + (k5 * N6) / 2 + (k6 * N9) / 2;
    Matrix[73] = (k5 * N14) / 2 - k4 * N13 + (k6 * N15) / 2;
    Matrix[74] = (k5 * N17) / 2 - k4 * N16 + (k6 * N18) / 2;
    Matrix[75] = k8 * N27 - k7 * N27;
    Matrix[76] = k8 * N24 - k7 * N23;
    Matrix[77] = k8 * N20 - k7 * N19;
}
//	Calculate element stress
void CBeam::ElementStress(double* stress, double* Displacement)
{
    const CBeamMaterial& material =
        static_cast<CBeamMaterial&>(*ElementMaterial); // Pointer to material of the element
    clear(stress, 3);
    double DX[3]; //	dx = x2-x1, dy = y2-y1, dz = z2-z1
    double L = 0; //	 beam length

    for (unsigned int i = 0; i < 3; i++)
    {
        DX[i] = nodes[1]->XYZ[i] - nodes[0]->XYZ[i];
    }

    L = sqrt(DX[0] * DX[0] + DX[1] * DX[1] + DX[2] * DX[2]);

    double S[6];
    for (unsigned int i = 0; i < 3; i++)
    {
        S[i] = -DX[i] * DX[i] * material.E / (L * L * L);
        S[i + 3] = -S[i];
    }

    for (unsigned int i = 0; i < 2; i++)
    {
        for (unsigned int j = 0; j < 3; j++)
        {
            if (LocationMatrix[i * 6 + j])
            {
                double a = S[i * 3 + j] * Displacement[LocationMatrix[i * 6 + j] - 1];
                stress[j] += S[i * 3 + j] * Displacement[LocationMatrix[i * 6 + j] - 1];
            }
        }
    }
}

#ifdef _VIB_
void CBeam::ElementMass(double* mass) {
}
#endif