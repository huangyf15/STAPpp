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

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//	Constructor
CInfEle::CInfEle()
{
    NEN = 4;	// Each element has 4 nodes
    nodes = new CNode*[NEN];
    
    ND = 12;
    LocationMatrix = new unsigned int[ND];

    ElementMaterial = nullptr;
}

//	Desconstructor
CInfEle::~CInfEle()
{
}

//	Read element data from stream Input
bool CInfEle::Read(ifstream& Input, unsigned int Ele, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int N;

    Input >> N;	// element number

    if (N != Ele + 1)
    {
        cerr << "*** Error *** Elements must be inputted in order !" << endl 
             << "   Expected element : " << Ele + 1 << endl
             << "   Provided element : " << N << endl;

        return false;
    }

    unsigned int MSet;	// Material property set number
    unsigned int N1, N2, N3, N4;	// numbers of 4 Nodes

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
	output << setw(5) << Ele+1 << setw(11) << nodes[0]->NodeNumber 
		   << setw(9) << nodes[1]->NodeNumber << setw(9) << nodes[2]->NodeNumber
           << setw(9) << nodes[3]->NodeNumber << setw(12) << ElementMaterial->nset << endl;
}

//  Generate location matrix: the global equation number that corresponding to each DOF of the element
//	Caution:  Equation number is numbered from 1 !
void CInfEle::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN; N++)
        for (unsigned int D = 0; D < 3; D++)
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

    const CInfiniteMaterial& material = static_cast<CInfiniteMaterial&>(*ElementMaterial);	// Pointer to material of the element

    double E = material.E;
    double neo = material.nu;  
    double eta, psi, Je;
    double J[2][2];
    double JT[2][2];
    double B[2][4];
    double Ke[36];
    double D[4];
    double x[4], y[4];
    for (unsigned i=0 ; i<36 ;i++){
        Ke[i]=0;
    }

    x[0] = nodes[0]->XYZ[0];
    y[0] = nodes[0]->XYZ[1];
    x[1] = nodes[1]->XYZ[0];
    y[1] = nodes[1]->XYZ[1];
    x[2] = 2 * x[0] - nodes[2]->XYZ[0];
    y[2] = 2 * y[0] - nodes[2]->XYZ[1];
    x[3] = 2 * x[1] - nodes[3]->XYZ[0];
    x[4] = 2 * x[1] - nodes[3]->XYZ[1];

    for (unsigned i=0 ; i<2 ;i++){
        for (unsigned j=0; j<2 ; j++){

            eta=pow(-1,i+1)/sqrt(3);
            psi=pow(-1,j+1)/sqrt(3); //Gauss point in local coordinate 

            J[0][0]=x[0]*(eta-1)/((1-psi)*(1-psi))+x[1]*(-1-eta)/((1-psi)*(1-psi))+x[2]*(1+eta)/((1-psi)*(1-psi))+x[3]*(1-eta)/((1-psi)*(1-psi));
            J[0][1]=y[0]*(eta-1)/((1-psi)*(1-psi))+y[1]*(-1-eta)/((1-psi)*(1-psi))+y[2]*(1+eta)/((1-psi)*(1-psi))+y[3]*(1-eta)/((1-psi)*(1-psi));
            J[1][0]=x[0]*psi/(1-psi)+x[1]*(-psi)/(1-psi)+x[2]*(1+psi)/(2*(1-psi))+x[3]*(-1-psi)/(2*(1-psi));
            J[1][1]=y[0]*psi/(1-psi)+y[1]*(-psi)/(1-psi)+y[2]*(1+psi)/(2*(1-psi))+y[3]*(-1-psi)/(2*(1-psi));
            Je=fabs(J[0][0]*J[1][1]-J[0][1]*J[1][0]);
            JT[0][0]=J[1][1]/Je;
            JT[0][1]=-J[0][1]/Je;
            JT[1][0]=-J[1][0]/Je;
            JT[1][1]=J[0][0]/Je;
            
            B[0][0]=JT[0][0]*(eta-1)/((1-psi)*(1-psi))+JT[0][1]*psi/(1-psi);
            B[0][1]=JT[0][0]*(-1-eta)/((1-psi)*(1-psi))+JT[0][1]*(-psi)/(1-psi);
            B[0][2]=JT[0][0]*(1+eta)/((1-psi)*(1-psi))+JT[0][1]*(psi+1)/(2*(1-psi));
            B[0][3]=JT[0][0]*(1-eta)/((1-psi)*(1-psi))+JT[0][1]*(-1-psi)/(2*(1-psi));
            B[1][0]=JT[1][0]*(eta-1)/((1-psi)*(1-psi))+JT[1][1]*psi/(1-psi);
            B[1][1]=JT[1][0]*(-1-eta)/((1-psi)*(1-psi))+JT[1][1]*(-psi)/(1-psi);
            B[1][2]=JT[1][0]*(1+eta)/((1-psi)*(1-psi))+JT[1][1]*(psi+1)/(2*(1-psi));
            B[1][3]=JT[1][0]*(1-eta)/((1-psi)*(1-psi))+JT[1][1]*(-1-psi)/(2*(1-psi));

            D[0]=E/(1-neo*neo);
            D[1]=neo*E/(1-neo*neo);
            D[2]=E/(1-neo*neo);
            D[3]=(1-neo)*E/(2-2*neo*neo);

            Ke[0]+=Je*(B[0][0]*B[0][0]*D[0]+B[1][0]*B[1][0]*D[3]);
            Ke[1]+=Je*(B[0][0]*B[1][0]*(D[1]+D[3]));
            Ke[2]+=Je*(B[1][0]*B[1][0]*D[2]+B[0][0]*B[0][0]*D[3]);
            Ke[3]+=Je*(B[0][0]*B[0][1]*D[0]+B[1][0]*B[1][1]*D[3]);
            Ke[4]+=Je*(B[1][0]*B[0][1]*D[1]+B[0][0]*B[1][1]*D[3]);
            Ke[5]+=Je*(B[0][1]*B[0][1]*D[0]+B[1][1]*B[1][1]*D[3]);
            Ke[6]+=Je*(B[0][0]*B[1][1]*D[1]+B[1][0]*B[0][1]*D[3]);
            Ke[7]+=Je*(B[1][0]*B[1][1]*D[2]+B[0][0]*B[0][1]*D[3]);
            Ke[8]+=Je*(B[0][1]*B[1][1]*D[1]+B[1][1]*B[0][1]*D[3]);
            Ke[9]+=Je*(B[1][1]*B[1][1]*D[2]+B[0][1]*B[0][1]*D[3]);
            Ke[10]+=Je*(B[0][0]*B[0][2]*D[0]+B[1][0]*B[1][2]*D[3]);
            Ke[11]+=Je*(B[1][0]*B[0][2]*D[1]+B[0][0]*B[1][2]*D[3]);
            Ke[12]+=Je*(B[0][1]*B[0][2]*D[0]+B[1][1]*B[1][2]*D[3]);
            Ke[13]+=Je*(B[1][1]*B[0][2]*D[1]+B[0][1]*B[1][2]*D[3]);
            Ke[14]+=Je*(B[0][2]*B[0][2]*D[0]+B[1][2]*B[1][2]*D[3]);
            Ke[15]+=Je*(B[0][0]*B[1][2]*D[1]+B[1][0]*B[0][2]*D[3]);
            Ke[16]+=Je*(B[1][0]*B[1][2]*D[2]+B[0][0]*B[0][2]*D[3]);
            Ke[17]+=Je*(B[0][1]*B[1][2]*D[1]+B[1][1]*B[0][2]*D[3]);
            Ke[18]+=Je*(B[1][1]*B[1][2]*D[2]+B[0][1]*B[0][2]*D[3]);
            Ke[19]+=Je*(B[0][2]*B[1][2]*D[1]+B[1][2]*B[0][2]*D[3]);
            Ke[20]+=Je*(B[1][2]*B[1][2]*D[2]+B[0][2]*B[0][2]*D[3]);
            Ke[21]+=Je*(B[0][0]*B[0][3]*D[0]+B[1][0]*B[1][3]*D[3]);
            Ke[22]+=Je*(B[1][0]*B[0][3]*D[1]+B[0][0]*B[1][3]*D[3]);
            Ke[23]+=Je*(B[0][1]*B[0][3]*D[0]+B[1][1]*B[1][3]*D[3]);
            Ke[24]+=Je*(B[1][1]*B[0][3]*D[1]+B[0][1]*B[1][3]*D[3]);
            Ke[25]+=Je*(B[0][2]*B[0][3]*D[0]+B[1][2]*B[1][3]*D[3]);
            Ke[26]+=Je*(B[1][2]*B[0][3]*D[1]+B[0][2]*B[1][3]*D[3]);
            Ke[27]+=Je*(B[0][3]*B[0][3]*D[0]+B[1][3]*B[1][3]*D[3]);
            Ke[28]+=Je*(B[0][0]*B[1][3]*D[1]+B[1][0]*B[0][3]*D[3]);
            Ke[29]+=Je*(B[1][0]*B[1][3]*D[2]+B[0][0]*B[0][3]*D[3]);
            Ke[30]+=Je*(B[0][1]*B[1][3]*D[1]+B[1][1]*B[0][3]*D[3]);
            Ke[31]+=Je*(B[1][1]*B[1][3]*D[2]+B[0][1]*B[0][3]*D[3]);
            Ke[32]+=Je*(B[0][2]*B[1][3]*D[1]+B[1][2]*B[0][3]*D[3]);
            Ke[33]+=Je*(B[1][2]*B[1][3]*D[2]+B[0][2]*B[0][3]*D[3]);
            Ke[34]+=Je*(B[0][3]*B[1][3]*D[1]+B[1][3]*B[0][3]*D[3]);
            Ke[35]+=Je*(B[1][3]*B[1][3]*D[2]+B[0][3]*B[0][3]*D[3]);
        }
    }
}

//	Calculate element stress 
void CInfEle::ElementStress(double* Q4stress, double* Displacement)
{
    clear(Q4stress,24);
    //	计算坐标变换
    double DX[9];		//	dx1 = x2-x1, dy1 = y2-y1, dz1 = z2-z1 , dx2=x3-x1 , dy2=y3-y1 , dz2=z3-z1
    for (unsigned int i = 0; i < 3; i++){
        DX[i] = nodes[1]->XYZ[i] - nodes[0]->XYZ[i];
        DX[i+3]=nodes[2]->XYZ[i] - nodes[0]->XYZ[i];
        DX[i+6]=nodes[3]->XYZ[i] - nodes[0]->XYZ[i];
    }
    
    double L[3],n[3][3],theta[2];//储存长度,4Q坐标系坐标轴在物理坐标系中的单位向量表示,以及31，41关于x'轴的夹角cos值（以21为x'轴）
    for (unsigned int i = 0; i < 3;i++){
        L[i]=sqrt(DX[3*i]*DX[3*i]+DX[3*i+1]*DX[3*i+1]+DX[3*i+2]*DX[3*i+2]);
        n[0][i]=DX[i]/L[0];
    }

    n[2][0]=DX[1]*DX[5]-DX[2]*DX[4];
    n[2][1]=DX[2]*DX[3]-DX[0]*DX[5];
    n[2][2]=DX[0]*DX[4]-DX[1]*DX[3];//叉乘得到z'轴向量
    double L2=sqrt(n[2][0]*n[2][0]+n[2][1]*n[2][1]+n[2][2]*n[2][2]);//归一化
    for (unsigned int i = 0; i<3 ;i++){
        n[2][i]=n[2][i]/L2;
    }

    n[1][0]=n[2][1]*n[0][2]-n[2][2]*n[0][1];
    n[1][1]=n[2][2]*n[0][0]-n[2][0]*n[0][2];
    n[1][2]=n[2][0]*n[0][1]-n[2][1]*n[0][0];//叉乘得到y'轴单位向量

    double x[4];
    double y[4];//储存4Q坐标系下节点的坐标
    x[0]=0 ; y[0]=0 ; x[1]=L[0] ; y[1]=0;
    x[2]=DX[3]*n[0][0]+DX[4]*n[0][1]+DX[5]*n[0][2];
    y[2]=DX[3]*n[1][0]+DX[4]*n[1][1]+DX[5]*n[1][2];
    x[3]=DX[6]*n[0][0]+DX[7]*n[0][1]+DX[8]*n[0][2];
    y[3]=DX[6]*n[1][0]+DX[7]*n[1][1]+DX[8]*n[1][2];
    
    CInfiniteMaterial* material = static_cast<CInfiniteMaterial*>(ElementMaterial);	// Pointer to material of the element

    double E = material->E;
    double nu=material->nu;  
    double eta,psi,Je;
    double J[2][2];
    double JT[2][2];
    double B[2][4];
    double Ke[36];
    double D[4];
    for (unsigned i=0; i<2 ;i++){
        for (unsigned j=0; j<2 ; j++){
            eta=pow(-1,i+1)/sqrt(3);
            psi=pow(-1,j+1)/sqrt(3);
            J[0][0]=x[0]*(eta-1)/4+x[1]*(1-eta)/4+x[2]*(1+eta)/4+x[3]*(-1-eta)/4;
            J[0][1]=y[0]*(eta-1)/4+y[1]*(1-eta)/4+y[2]*(1+eta)/4+y[3]*(-1-eta)/4;
            J[1][0]=x[0]*(psi-1)/4+x[1]*(-1-psi)/4+x[2]*(1+psi)/4+x[3]*(1-psi)/4;
            J[1][1]=y[0]*(psi-1)/4+y[1]*(-1-psi)/4+y[2]*(1+psi)/4+y[3]*(1-psi)/4;
            Je=J[0][0]*J[1][1]-J[0][1]*J[1][0];
            JT[0][0]=J[1][1]/Je;
            JT[0][1]=-J[0][1]/Je;
            JT[1][0]=-J[1][0]/Je;
            JT[1][1]=J[0][0]/Je;
            
            B[0][0]=JT[0][0]*(eta-1)/4+JT[0][1]*(psi-1)/4;
            B[0][1]=JT[0][0]*(1-eta)/4+JT[0][1]*(-psi-1)/4;
            B[0][2]=JT[0][0]*(eta+1)/4+JT[0][1]*(psi+1)/4;
            B[0][3]=JT[0][0]*(-eta-1)/4+JT[0][1]*(1-psi)/4;
            B[1][0]=JT[1][0]*(eta-1)/4+JT[1][1]*(psi-1)/4;
            B[1][1]=JT[1][0]*(1-eta)/4+JT[1][1]*(-psi-1)/4;
            B[1][2]=JT[1][0]*(eta+1)/4+JT[1][1]*(psi+1)/4;
            B[1][3]=JT[1][0]*(-eta-1)/4+JT[1][1]*(1-psi)/4;

            D[0]=E/(1-nu*nu);
            D[1]=nu*E/(1-nu*nu);
            D[2]=E/(1-nu*nu);
            D[3]=(1-nu)*E/(2-2*nu);

            double de[8];//单元位移
            for (unsigned k=0; k<8 ; k++){
                de[k]=0;
            }//清零
            for (unsigned k=0; k<4 ; k++){
                for (unsigned v=0 ; v<3 ; v++){
                    if (LocationMatrix[3*k+v]){
                        de[2*k]  += n[0][v]*Displacement[LocationMatrix[3*k+v]-1];
                        de[2*k+1]+= n[1][v]*Displacement[LocationMatrix[3*k+v]-1];
                    }
                }  
            }

            //计算高斯点的应力
            double strain[3],stress[3];//应变，平面应力作为中间变量
            strain[0]=B[0][0]*de[0]+B[0][1]*de[2]+B[0][2]*de[4]+B[0][3]*de[6];
            strain[1]=B[1][0]*de[1]+B[1][1]*de[3]+B[1][2]*de[5]+B[1][3]*de[7];
            strain[2]=B[1][0]*de[0]+B[0][0]*de[1]+B[1][1]*de[2]+B[0][1]*de[3]+B[1][2]*de[4]+B[0][2]*de[5]+B[1][3]*de[6]+B[0][3]*de[7];
            stress[0]=D[0]*strain[0]+D[1]*strain[1];
            stress[1]=D[1]*strain[0]+D[2]*strain[1];
            stress[2]=D[3]*strain[2];
            if ((!i)&&(!j)){
                Q4stress[0]=n[0][0]*n[0][0]*stress[0]+2*n[0][0]*n[1][0]*stress[2]+n[1][0]*n[1][0]*stress[1];
                Q4stress[1]=n[0][1]*n[0][1]*stress[0]+2*n[0][1]*n[1][1]*stress[2]+n[1][1]*n[1][1]*stress[1];
                Q4stress[2]=n[0][2]*n[0][2]*stress[0]+2*n[0][2]*n[1][2]*stress[2]+n[1][2]*n[1][2]*stress[1];
                Q4stress[3]=n[0][1]*n[0][2]*stress[0]+(n[1][1]*n[0][2]+n[0][1]*n[1][2])*stress[2]+n[1][1]*n[1][2]*stress[1];
                Q4stress[4]=n[0][0]*n[0][2]*stress[0]+(n[1][0]*n[0][2]+n[0][0]*n[1][2])*stress[2]+n[1][0]*n[1][2]*stress[1];
                Q4stress[5]=n[0][0]*n[0][1]*stress[0]+(n[1][0]*n[0][1]+n[0][0]*n[1][1])*stress[2]+n[1][0]*n[1][1]*stress[1];
            }

            if ((!i)&&(j)){
                Q4stress[6]=n[0][0]*n[0][0]*stress[0]+2*n[0][0]*n[1][0]*stress[2]+n[1][0]*n[1][0]*stress[1];
                Q4stress[7]=n[0][1]*n[0][1]*stress[0]+2*n[0][1]*n[1][1]*stress[2]+n[1][1]*n[1][1]*stress[1];
                Q4stress[8]=n[0][2]*n[0][2]*stress[0]+2*n[0][2]*n[1][2]*stress[2]+n[1][2]*n[1][2]*stress[1];
                Q4stress[9]=n[0][1]*n[0][2]*stress[0]+(n[1][1]*n[0][2]+n[0][1]*n[1][2])*stress[2]+n[1][1]*n[1][2]*stress[1];
                Q4stress[10]=n[0][0]*n[0][2]*stress[0]+(n[1][0]*n[0][2]+n[0][0]*n[1][2])*stress[2]+n[1][0]*n[1][2]*stress[1];
                Q4stress[11]=n[0][0]*n[0][1]*stress[0]+(n[1][0]*n[0][1]+n[0][0]*n[1][1])*stress[2]+n[1][0]*n[1][1]*stress[1];
            }        

            if ((i)&&(!j)){
                Q4stress[12]=n[0][0]*n[0][0]*stress[0]+2*n[0][0]*n[1][0]*stress[2]+n[1][0]*n[1][0]*stress[1];
                Q4stress[13]=n[0][1]*n[0][1]*stress[0]+2*n[0][1]*n[1][1]*stress[2]+n[1][1]*n[1][1]*stress[1];
                Q4stress[14]=n[0][2]*n[0][2]*stress[0]+2*n[0][2]*n[1][2]*stress[2]+n[1][2]*n[1][2]*stress[1];
                Q4stress[15]=n[0][1]*n[0][2]*stress[0]+(n[1][1]*n[0][2]+n[0][1]*n[1][2])*stress[2]+n[1][1]*n[1][2]*stress[1];
                Q4stress[16]=n[0][0]*n[0][2]*stress[0]+(n[1][0]*n[0][2]+n[0][0]*n[1][2])*stress[2]+n[1][0]*n[1][2]*stress[1];
                Q4stress[17]=n[0][0]*n[0][1]*stress[0]+(n[1][0]*n[0][1]+n[0][0]*n[1][1])*stress[2]+n[1][0]*n[1][1]*stress[1];
            }  

            if ((i)&&(j)){
                Q4stress[18]=n[0][0]*n[0][0]*stress[0]+2*n[0][0]*n[1][0]*stress[2]+n[1][0]*n[1][0]*stress[1];
                Q4stress[19]=n[0][1]*n[0][1]*stress[0]+2*n[0][1]*n[1][1]*stress[2]+n[1][1]*n[1][1]*stress[1];
                Q4stress[20]=n[0][2]*n[0][2]*stress[0]+2*n[0][2]*n[1][2]*stress[2]+n[1][2]*n[1][2]*stress[1];
                Q4stress[21]=n[0][1]*n[0][2]*stress[0]+(n[1][1]*n[0][2]+n[0][1]*n[1][2])*stress[2]+n[1][1]*n[1][2]*stress[1];
                Q4stress[22]=n[0][0]*n[0][2]*stress[0]+(n[1][0]*n[0][2]+n[0][0]*n[1][2])*stress[2]+n[1][0]*n[1][2]*stress[1];
                Q4stress[23]=n[0][0]*n[0][1]*stress[0]+(n[1][0]*n[0][1]+n[0][0]*n[1][1])*stress[2]+n[1][0]*n[1][1]*stress[1];
            }  
        }
    }
}
void CInfEle::ElementGauss(double* Coordinate){
    clear(Coordinate,12);
    //	计算坐标变换
    double DX[12];		//	dx1 = x2-x1, dy1 = y2-y1, dz1 = z2-z1 , dx2=x3-x1 , dy2=y3-y1 , dz2=z3-z1
    for (unsigned int i = 0; i < 3; i++){
        DX[i] = nodes[1]->XYZ[i] - nodes[0]->XYZ[i];
        DX[i+3]=nodes[2]->XYZ[i] - nodes[0]->XYZ[i];
        DX[i+6]=nodes[3]->XYZ[i] - nodes[0]->XYZ[i];
        DX[i+9]=nodes[0]->XYZ[i];
    }
    
    double L[3],n[3][3],theta[2];//储存长度,4Q坐标系坐标轴在物理坐标系中的单位向量表示,以及31，41关于x'轴的夹角cos值（以21为x'轴）
    for (unsigned int i = 0; i < 3;i++){
        L[i]=sqrt(DX[3*i]*DX[3*i]+DX[3*i+1]*DX[3*i+1]+DX[3*i+2]*DX[3*i+2]);
        n[0][i]=DX[i]/L[0];
    }

    n[2][0]=DX[1]*DX[5]-DX[2]*DX[4];
    n[2][1]=DX[2]*DX[3]-DX[0]*DX[5];
    n[2][2]=DX[0]*DX[4]-DX[1]*DX[3];//叉乘得到z'轴向量
    double L2=sqrt(n[2][0]*n[2][0]+n[2][1]*n[2][1]+n[2][2]*n[2][2]);//归一化
    for (unsigned int i = 0; i<3 ;i++){
        n[2][i]=n[2][i]/L2;
    }

    n[1][0]=n[2][1]*n[0][2]-n[2][2]*n[0][1];
    n[1][1]=n[2][2]*n[0][0]-n[2][0]*n[0][2];
    n[1][2]=n[2][0]*n[0][1]-n[2][1]*n[0][0];//叉乘得到y'轴单位向量

    double x[4];
    double y[4];//储存4Q坐标系下节点的坐标
    x[0]=0 ; y[0]=0 ; x[1]=L[0] ; y[1]=0;
    x[2]=DX[3]*n[0][0]+DX[4]*n[0][1]+DX[5]*n[0][2];
    y[2]=DX[3]*n[1][0]+DX[4]*n[1][1]+DX[5]*n[1][2];
    x[3]=DX[6]*n[0][0]+DX[7]*n[0][1]+DX[8]*n[0][2];
    y[3]=DX[6]*n[1][0]+DX[7]*n[1][1]+DX[8]*n[1][2];

    double eta,psi,C1,C2;
    for (unsigned i=0; i<2 ;i++){
        for (unsigned j=0; j<2 ; j++){
            eta=pow(-1,i+1)/sqrt(3);
            psi=pow(-1,j+1)/sqrt(3); 
            if ((!i)&(!j)){
                C1=(1-psi)*(1-eta)*x[0]/4+(1+psi)*(1-eta)*x[1]/4+(1+psi)*(1+eta)*x[2]/4+(1-psi)*(1+eta)*x[3]/4;
                C2=(1-psi)*(1-eta)*y[0]/4+(1+psi)*(1-eta)*y[1]/4+(1+psi)*(1+eta)*y[2]/4+(1-psi)*(1+eta)*y[3]/4;
                Coordinate[0]=n[0][0]*C1+n[1][0]*C2+DX[9];
                Coordinate[1]=n[0][1]*C1+n[1][1]*C2+DX[10];
                Coordinate[2]=n[0][2]*C1+n[1][2]*C2+DX[11];
            }

            if ((!i)&(j)){
                C1=(1-psi)*(1-eta)*x[0]/4+(1+psi)*(1-eta)*x[1]/4+(1+psi)*(1+eta)*x[2]/4+(1-psi)*(1+eta)*x[3]/4;
                C2=(1-psi)*(1-eta)*y[0]/4+(1+psi)*(1-eta)*y[1]/4+(1+psi)*(1+eta)*y[2]/4+(1-psi)*(1+eta)*y[3]/4;
                Coordinate[3]=n[0][0]*C1+n[1][0]*C2+DX[9];
                Coordinate[4]=n[0][1]*C1+n[1][1]*C2+DX[10];
                Coordinate[5]=n[0][2]*C1+n[1][2]*C2+DX[11];
            }

            if ((i)&(!j)){
                C1=(1-psi)*(1-eta)*x[0]/4+(1+psi)*(1-eta)*x[1]/4+(1+psi)*(1+eta)*x[2]/4+(1-psi)*(1+eta)*x[3]/4;
                C2=(1-psi)*(1-eta)*y[0]/4+(1+psi)*(1-eta)*y[1]/4+(1+psi)*(1+eta)*y[2]/4+(1-psi)*(1+eta)*y[3]/4;
                Coordinate[6]=n[0][0]*C1+n[1][0]*C2+DX[9];
                Coordinate[7]=n[0][1]*C1+n[1][1]*C2+DX[10];
                Coordinate[8]=n[0][2]*C1+n[1][2]*C2+DX[11];
            }

            if ((i)&(j)){
                C1=(1-psi)*(1-eta)*x[0]/4+(1+psi)*(1-eta)*x[1]/4+(1+psi)*(1+eta)*x[2]/4+(1-psi)*(1+eta)*x[3]/4;
                C2=(1-psi)*(1-eta)*y[0]/4+(1+psi)*(1-eta)*y[1]/4+(1+psi)*(1+eta)*y[2]/4+(1-psi)*(1+eta)*y[3]/4;
                Coordinate[9]=n[0][0]*C1+n[1][0]*C2+DX[9];
                Coordinate[10]=n[0][1]*C1+n[1][1]*C2+DX[10];
                Coordinate[11]=n[0][2]*C1+n[1][2]*C2+DX[11];
            }
        }
    } 
}

