/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.02, October 27, 2017                                        */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Elements/Infinite_4Q.h"

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace std;

//	Constructor
CInfEle::CInfEle()
{
    NEN = 4; // Each element has 4 nodes
    nodes = new CNode*[NEN];

    ND = 8;
    LocationMatrix = new unsigned int[ND];

    ElementMaterial = nullptr;
}

//	Desconstructor
CInfEle::~CInfEle() {}

//	Read element data from stream Input
bool CInfEle::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int N;

    Input >> N; // element number

    if (N != Ele + 1)
    {
        cerr << "*** Error *** Elements must be inputted in order !" << endl
             << "   Expected element : " << Ele + 1 << endl
             << "   Provided element : " << N << endl;

        return false;
    }

    unsigned int MSet;           // Material property set number
    unsigned int N1, N2, N3, N4; // numbers of 4 Nodes

    Input >> N1 >> N2 >> N3 >> N4 >> MSet;
    ElementMaterial = static_cast<CInfiniteMaterial*>(MaterialSets) + MSet - 1;
    nodes[0] = &NodeList[N1 - 1];
    nodes[1] = &NodeList[N2 - 1];
    nodes[2] = &NodeList[N3 - 1];
    nodes[3] = &NodeList[N4 - 1];

    return true;
}

//	Write element data to stream OutputFile
void CInfEle::Write(COutputter& output, unsigned int Ele)
{
    output << setw(5) << Ele + 1 << setw(11) << nodes[0]->NodeNumber << setw(9)
           << nodes[1]->NodeNumber << setw(9) << nodes[2]->NodeNumber << setw(9)
           << nodes[3]->NodeNumber << setw(12) << ElementMaterial->nset << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the
//  element
//	Caution:  Equation number is numbered from 1 !
void CInfEle::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN; N++)
        for (unsigned int D = 0; D < 2; D++)
            LocationMatrix[i++] = nodes[N]->bcode[D];
}

//	Return the size of the element stiffness matrix (stored as an array column by column)
//	For infinite 4Q element, element stiffness is a 8x8 matrix, whose upper triangular part
//	has 36 elements
unsigned int CInfEle::SizeOfStiffnessMatrix() { return 36; }

//	Calculate element stiffness matrix
//	Upper triangular matrix, stored as an array column by colum starting from the diagonal element
void CInfEle::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());

    //	Calculate element stiffness matrix

    const CInfiniteMaterial& material =
        static_cast<CInfiniteMaterial&>(*ElementMaterial); // Pointer to material of the element

    double E = material.E;
    double nu = material.nu;
    double eta, psi, Je;
    double J[2][2];
    double JT[2][2];
    double B[2][4];
    double Ke[36];
    double D[4];
    double x[4], y[4];
    for (unsigned i = 0; i < 36; i++)
    {
        Ke[i] = 0;
    }

    x[0] = nodes[0]->XYZ[0];
    y[0] = nodes[0]->XYZ[1];
    x[1] = nodes[1]->XYZ[0];
    y[1] = nodes[1]->XYZ[1];
    x[2] = nodes[2]->XYZ[0];
    y[2] = nodes[2]->XYZ[1];
    x[3] = nodes[3]->XYZ[0];
    y[3] = nodes[3]->XYZ[1];

    for (unsigned i = 0; i < 2; i++)
    {
        for (unsigned j = 0; j < 2; j++)
        {

            eta = pow(-1, i + 1) / sqrt(3);
            psi = pow(-1, j + 1) / sqrt(3); // Gauss point in local coordinate

            J[0][0] = x[0] * (eta - 1) / ((1 - psi) * (1 - psi)) +
                      x[1] * (-1 - eta) / ((1 - psi) * (1 - psi)) +
                      x[2] * (1 + eta) / ((1 - psi) * (1 - psi)) +
                      x[3] * (1 - eta) / ((1 - psi) * (1 - psi));
            J[0][1] = y[0] * (eta - 1) / ((1 - psi) * (1 - psi)) +
                      y[1] * (-1 - eta) / ((1 - psi) * (1 - psi)) +
                      y[2] * (1 + eta) / ((1 - psi) * (1 - psi)) +
                      y[3] * (1 - eta) / ((1 - psi) * (1 - psi));
            J[1][0] = x[0] * psi / (1 - psi) + x[1] * (-psi) / (1 - psi) +
                      x[2] * (1 + psi) / (2 * (1 - psi)) + x[3] * (-1 - psi) / (2 * (1 - psi));
            J[1][1] = y[0] * psi / (1 - psi) + y[1] * (-psi) / (1 - psi) +
                      y[2] * (1 + psi) / (2 * (1 - psi)) + y[3] * (-1 - psi) / (2 * (1 - psi));
            Je = fabs(J[0][0] * J[1][1] - J[0][1] * J[1][0]);
            JT[0][0] = J[1][1] / Je;
            JT[0][1] = -J[0][1] / Je;
            JT[1][0] = -J[1][0] / Je;
            JT[1][1] = J[0][0] / Je;

            B[0][0] = JT[0][0] * (eta - 1) / ((1 - psi) * (1 - psi)) + JT[0][1] * psi / (1 - psi);
            B[0][1] =
                JT[0][0] * (-1 - eta) / ((1 - psi) * (1 - psi)) + JT[0][1] * (-psi) / (1 - psi);
            B[0][2] = JT[0][0] * (1 + eta) / ((1 - psi) * (1 - psi)) +
                      JT[0][1] * (psi + 1) / (2 * (1 - psi));
            B[0][3] = JT[0][0] * (1 - eta) / ((1 - psi) * (1 - psi)) +
                      JT[0][1] * (-1 - psi) / (2 * (1 - psi));
            B[1][0] = JT[1][0] * (eta - 1) / ((1 - psi) * (1 - psi)) + JT[1][1] * psi / (1 - psi);
            B[1][1] =
                JT[1][0] * (-1 - eta) / ((1 - psi) * (1 - psi)) + JT[1][1] * (-psi) / (1 - psi);
            B[1][2] = JT[1][0] * (1 + eta) / ((1 - psi) * (1 - psi)) +
                      JT[1][1] * (psi + 1) / (2 * (1 - psi));
            B[1][3] = JT[1][0] * (1 - eta) / ((1 - psi) * (1 - psi)) +
                      JT[1][1] * (-1 - psi) / (2 * (1 - psi));

            D[0] = E / (1 - nu * nu);
            D[1] = nu * E / (1 - nu * nu);
            D[2] = E / (1 - nu * nu);
            D[3] = (1 - nu) * E / (2 - 2 * nu * nu);

            Matrix[0] += Je * (B[0][0] * B[0][0] * D[0] + B[1][0] * B[1][0] * D[3]);
            Matrix[1] += Je * (B[1][0] * B[1][0] * D[2] + B[0][0] * B[0][0] * D[3]);
            Matrix[2] += Je * (B[0][0] * B[1][0] * (D[1] + D[3]));
            Matrix[3] += Je * (B[0][1] * B[0][1] * D[0] + B[1][1] * B[1][1] * D[3]);
            Matrix[4] += Je * (B[1][0] * B[0][1] * D[1] + B[0][0] * B[1][1] * D[3]);
            Matrix[5] += Je * (B[0][0] * B[0][1] * D[0] + B[1][0] * B[1][1] * D[3]);
            Matrix[6] += Je * (B[1][1] * B[1][1] * D[2] + B[0][1] * B[0][1] * D[3]);
            Matrix[7] += Je * (B[0][1] * B[1][1] * D[1] + B[1][1] * B[0][1] * D[3]);
            Matrix[8] += Je * (B[1][0] * B[1][1] * D[2] + B[0][0] * B[0][1] * D[3]);
            Matrix[9] += Je * (B[0][0] * B[1][1] * D[1] + B[1][0] * B[0][1] * D[3]);
            Matrix[10] += Je * (B[0][2] * B[0][2] * D[0] + B[1][2] * B[1][2] * D[3]);
            Matrix[11] += Je * (B[1][1] * B[0][2] * D[1] + B[0][1] * B[1][2] * D[3]);
            Matrix[12] += Je * (B[0][1] * B[0][2] * D[0] + B[1][1] * B[1][2] * D[3]);
            Matrix[13] += Je * (B[1][0] * B[0][2] * D[1] + B[0][0] * B[1][2] * D[3]);
            Matrix[14] += Je * (B[0][0] * B[0][2] * D[0] + B[1][0] * B[1][2] * D[3]);
            Matrix[15] += Je * (B[1][2] * B[1][2] * D[2] + B[0][2] * B[0][2] * D[3]);
            Matrix[16] += Je * (B[0][2] * B[1][2] * D[1] + B[1][2] * B[0][2] * D[3]);
            Matrix[17] += Je * (B[1][1] * B[1][2] * D[2] + B[0][1] * B[0][2] * D[3]);
            Matrix[18] += Je * (B[0][1] * B[1][2] * D[1] + B[1][1] * B[0][2] * D[3]);
            Matrix[19] += Je * (B[1][0] * B[1][2] * D[2] + B[0][0] * B[0][2] * D[3]);
            Matrix[20] += Je * (B[0][0] * B[1][2] * D[1] + B[1][0] * B[0][2] * D[3]);
            Matrix[21] += Je * (B[0][3] * B[0][3] * D[0] + B[1][3] * B[1][3] * D[3]);
            Matrix[22] += Je * (B[1][2] * B[0][3] * D[1] + B[0][2] * B[1][3] * D[3]);
            Matrix[23] += Je * (B[0][2] * B[0][3] * D[0] + B[1][2] * B[1][3] * D[3]);
            Matrix[24] += Je * (B[1][1] * B[0][3] * D[1] + B[0][1] * B[1][3] * D[3]);
            Matrix[25] += Je * (B[0][1] * B[0][3] * D[0] + B[1][1] * B[1][3] * D[3]);
            Matrix[26] += Je * (B[1][0] * B[0][3] * D[1] + B[0][0] * B[1][3] * D[3]);
            Matrix[27] += Je * (B[0][0] * B[0][3] * D[0] + B[1][0] * B[1][3] * D[3]);
            Matrix[28] += Je * (B[1][3] * B[1][3] * D[2] + B[0][3] * B[0][3] * D[3]);
            Matrix[29] += Je * (B[0][3] * B[1][3] * D[1] + B[1][3] * B[0][3] * D[3]);
            Matrix[30] += Je * (B[1][2] * B[1][3] * D[2] + B[0][2] * B[0][3] * D[3]);
            Matrix[31] += Je * (B[0][2] * B[1][3] * D[1] + B[1][2] * B[0][3] * D[3]);
            Matrix[32] += Je * (B[1][1] * B[1][3] * D[2] + B[0][1] * B[0][3] * D[3]);
            Matrix[33] += Je * (B[0][1] * B[1][3] * D[1] + B[1][1] * B[0][3] * D[3]);
            Matrix[34] += Je * (B[1][0] * B[1][3] * D[2] + B[0][0] * B[0][3] * D[3]);
            Matrix[35] += Je * (B[0][0] * B[1][3] * D[1] + B[1][0] * B[0][3] * D[3]);
        }
    }
}

//	Calculate element stress
void CInfEle::ElementStress(double* Q4stress, double* Displacement) {}
void CInfEle::ElementGauss(double* Coordinate) {}
