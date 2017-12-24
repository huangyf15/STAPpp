#include "Elements/8H.h"
#include "Eigen/Dense"
#include <iostream>
#include <iomanip>
#include <cmath>
#include "SPR8H.h"

using namespace Eigen;
using namespace std;
typedef short unsigned int s_uint;
typedef unsigned int uint;

void  CHex::ElementPostSPR(double* stressG, double* Displacement , double* PrePositions, double* PostPositions, double* PositionG)
{
	

//position before displacement for nodes
for (unsigned int i=0;i<8;i++)
	{ 
		for (unsigned int j=0;j<3;j++)
	   {
			PrePositions[i*3+j]=nodes[i]->XYZ[j];
	    }
	}

	

//ideal displacement matrix de
	VectorXd idealdisp(24); 
		for (unsigned int i=0;i<24;i++ )
		{ 
			if (LocationMatrix[i])
				{idealdisp(i) = Displacement[LocationMatrix[i]-1];}
			else
				{idealdisp(i) = 0.0;}

		
			PostPositions[i] = PrePositions[i] + idealdisp(i);

		}
	
	CHexMaterial* material = dynamic_cast<CHexMaterial*>(ElementMaterial);	// Pointer to material of the element
	
	//construct De
	double v = material->nu;
	double k = material->E * (1-v)/(1+v)/(1-2*v);

	MatrixXd D(6,6); // constitutive matrix
	D << 1,			v/(1-v),	v/(1-v),		0,					0,					0,
		v/(1-v),	1,			v/(1-v),		0,					0,					0,
		v/(1-v),	v/(1-v),    1,				0,					0,					0,
		0,			0,		    0,				(1-2*v)/2/(1-v),    0,					0,
		0,			0,			0,				0,					(1-2*v)/2/(1-v),	0,
		0,			0,			0,				0,					0,					(1-2*v)/2/(1-v);
	D=D*k;

	
	//construct coordinate matrix
	MatrixXd coorxyz(8,3);
	for (unsigned int i=0;i<8;i++)
	{ 
		for (unsigned int j=0;j<3;j++)
	   {
			coorxyz(i,j)=nodes[i]->XYZ[j];
	    }
	}


// construct jacobi matrix
	double xi;
	double eta;
	double zeta;
	double detj;

	double xi8[8];
	double eta8[8];
	double zeta8[8];

	xi8[0]=0.577350269189626;
	xi8[1]=0.577350269189626;
	xi8[2]=-0.577350269189626;
	xi8[3]=-0.577350269189626;
	xi8[4]=0.577350269189626;
	xi8[5]=0.577350269189626;
	xi8[6]=-0.577350269189626;
	xi8[7]=-0.577350269189626;
	
	eta8[0]=-0.577350269189626;
	eta8[1]=0.577350269189626;
	eta8[2]=0.577350269189626;
	eta8[3]=-0.577350269189626;
	eta8[4]=-0.577350269189626;
	eta8[5]=0.577350269189626;
	eta8[6]=0.577350269189626;
	eta8[7]=-0.577350269189626;

	zeta8[0]=-0.577350269189626;
	zeta8[1]=-0.577350269189626;
	zeta8[2]=-0.577350269189626;
	zeta8[3]=-0.577350269189626;
	zeta8[4]=0.577350269189626;
	zeta8[5]=0.577350269189626;
	zeta8[6]=0.577350269189626;
	zeta8[7]=0.577350269189626;

	
	MatrixXd GN(3,8);
	MatrixXd J;
	MatrixXd Jni;
	MatrixXd Bele; // elements in Be
	MatrixXd Be(6,24);
	MatrixXd stressXYZ(6,8); //8 gauss point, 6 stress 
	// shape function 
	MatrixXd N(1,8);
	// coordinate for gauss points
	MatrixXd coorG;
	
	for (unsigned p=0;p<8;p++)
	{
		
		xi   = xi8[p];
		eta  = eta8[p];
		zeta = zeta8[p];
		

		GN << (1-eta)*(1-zeta), (1+eta)*(1-zeta),  -(1+eta)*(1-zeta),  -(1-eta)*(1-zeta),  (1-eta)*(1+zeta),  (1+eta)*(1+zeta),  -(1+eta)*(1+zeta),  -(1-eta)*(1+zeta),  
			  -(1+xi)*(1-zeta),  (1+xi)*(1-zeta),    (1-xi)*(1-zeta),  -(1-xi)*(1-zeta),   -(1+xi)*(1+zeta),   (1+xi)*(1+zeta),   (1-xi)*(1+zeta),   -(1-xi)*(1+zeta),  
			  -(1+xi)*(1-eta),   -(1+xi)*(1+eta),    -(1-xi)*(1+eta),   -(1-xi)*(1-eta),    (1+xi)*(1-eta),    (1+xi)*(1+eta),    (1-xi)*(1+eta),    (1-xi)*(1-eta);    
		GN=GN/8; // coefficient

		J=GN*coorxyz;
		Jni=J.inverse();
		detj=J.determinant();
		Bele=Jni*GN;
	

		// assign the value of Be	
		Be << Bele(0,0),0,0,				Bele(0,1), 0,0,				Bele(0,2), 0,0,				Bele(0,3), 0,0,				Bele(0,4), 0,0,				Bele(0,5), 0,0,				Bele(0,6), 0,0,				Bele(0,7), 0,0,
			  0,Bele(1,0),0,				0,Bele(1,1),0,				0,Bele(1,2),0,				0,Bele(1,3),0,				0,Bele(1,4),0,				0,Bele(1,5),0,				0,Bele(1,6),0,				0,Bele(1,7),0,
			  0,0,Bele(2,0),				0,0,Bele(2,1),				0,0,Bele(2,2),				0,0,Bele(2,3),				0,0,Bele(2,4),				0,0,Bele(2,5),				0,0,Bele(2,6),				0,0,Bele(2,7),
			  Bele(1,0),Bele(0,0),0,		Bele(1,1),Bele(0,1),0,		Bele(1,2),Bele(0,2),0,		Bele(1,3),Bele(0,3),0,		Bele(1,4),Bele(0,4),0,		Bele(1,5),Bele(0,5),0,		Bele(1,6),Bele(0,6),0,		Bele(1,7),Bele(0,7),0,
			  0,Bele(2,0),Bele(1,0),		0,Bele(2,1),Bele(1,1),		0,Bele(2,2),Bele(1,2),		0,Bele(2,3),Bele(1,3),		0,Bele(2,4),Bele(1,4),		0,Bele(2,5),Bele(1,5),		0,Bele(2,6),Bele(1,6),		0,Bele(2,7),Bele(1,7),
			  Bele(2,0),0,Bele(0,0),		Bele(2,1),0,Bele(0,1),		Bele(2,2),0,Bele(0,2),		Bele(2,3),0,Bele(0,3),		Bele(2,4),0,Bele(0,4),		Bele(2,5),0,Bele(0,5),		Bele(2,6),0,Bele(0,6),		Bele(2,7),0,Bele(0,7);
	

		// 6,1 for each gauss point stress 
		//a column of stressXYZ
		stressXYZ.col(p)=D*Be*idealdisp;  
		//loop for every gauss point 

		// stress for gauss points
		for (unsigned int i=0;i<6;i++)
		{
			stressG[p*6+i]=stressXYZ(i,p);
		}
	}
	
// calculate the coordinate of gauss points

	for (unsigned int i=0; i<8; i++)
	{
		xi   = xi8[i];
		eta  = eta8[i];
		zeta = zeta8[i];
		N << (1+xi)*(1-eta)*(1-zeta), (1+xi)*(1+eta)*(1-zeta),  (1-xi)*(1+eta)*(1-zeta),  (1-xi)*(1-eta)*(1-zeta),  (1+xi)*(1-eta)*(1+zeta),  (1+xi)*(1+eta)*(1+zeta),  (1-xi)*(1+eta)*(1+zeta),  (1-xi)*(1-eta)*(1+zeta);
		N=N/8; // coefficient

		coorG=N*coorxyz;

		for(unsigned int j=0;j<3;j++)
		{
			PositionG[i*3+j]=coorG(j);

		}
	}
}

void StressSPR(double* stress_SPR, double* stressG, double* PrePositions, double* PositionG, 
			   uint* Ele_NodeNumber, uint NUME, uint NUMNP)
{
	uint* NNE = new uint[NUMNP]; // number of neibourhood elements
	s_uint* CNT = new s_uint[NUMNP]; // how many times the stress of the node has been calculated
	double* stress_tmp = new double [NUMNP*6]; //按全局节点号存储的应力值
	// 先清点每个节点连接的单元个数
	for(uint Np = 0; Np < NUMNP; Np++)
	{
		CNT[Np] = 0;
		NNE[Np] = 0;
		for(uint Ele = 0; Ele < NUME; Ele++)
		{
			for(uint N = 0; N < 8; N++)
			{
				if(Ele_NodeNumber[Ele*8+N] == Np+1)
					NNE[Np]++;
				else
					continue;
			}
		}
	}

	// 按节点遍历，如果与该节点相连的单元个数大于等于2（NNE>=2），则将这些单元作为SPR的patch。
	// 除计算该节点处的应力外，还计算该patch中NNE=1的节点的应力
	// 最后NNE=1的节点的应力取按各个patch计算结果的平均
	for (uint i=1; i<NUMNP*6; i++)
		stress_tmp[i] = 0;    //对应力进行初始化

	for(uint Np = 0; Np < NUMNP; Np++)
	{
		if(NNE[Np]>=2)
		{
			uint* Ele_Np = new uint[NNE[Np]];     // 与Np节点相连的单元的编号
			s_uint index = 0; // 与Np节点相连的单元计数
			for(uint Ele = 0; Ele < NUME; Ele++ ) // 找出与Np节点相连的单元编号	
			{
				for(s_uint N = 0; N < 8; N++)
				{
					if(Ele_NodeNumber[N+Ele*8] == Np+1)
					{
						Ele_Np[index] = Ele;
						index++;
						break;
					}
				}
				if(index>=NNE[Np]) 
					break;
			}

			s_uint N_SPR = 1; //patch中需要计算应力的点的个数，Np已经算一个，NNE等于1的点也要算
			uint* NP_SPR = new uint [NNE[Np]*8]; 
			//patch中需要计算应力的节点的编号，因为事先不知道有多少个点要计算的，所以大小取了富余
			NP_SPR[0] = Np+1; //Np是patch的中心点，是第一个需要计算应力的点；
			for(s_uint i=0;i<NNE[Np];i++)
			{
				for(s_uint N=0;N<8;N++)
				{
					if (NNE[Ele_NodeNumber[Ele_Np[i]*8+N]]==1)
					{
						NP_SPR[N_SPR] = Ele_NodeNumber[Ele_Np[i]*8+N];
						N_SPR++;
					}					
				}

			}

            for (s_uint k=0;k<N_SPR;k++)
			{
				CNT[NP_SPR[k]-1]++;  //节点应力每被计算一次，计数增加1
			}

			//用二阶完备多项式拟合应力场，基函数为1,x,y,z,x^2,y^2,z^2,xy,xz,yz
		
			for (s_uint j=0;j<6;j++) //按6个应力分量遍历
			{
				// 先提取该patch中计算所需要的数据
				// a = inv(A)*b A=sum(ppT)
				MatrixXd A(10,10); A=MatrixXd::Zero(10,10);
				MatrixXd b(10,1);  b=MatrixXd::Zero(10,1);
				MatrixXd p(10,1);  
				MatrixXd a(10,1);  //二阶完备多项式的系数
				for (s_uint k=0;k<NNE[Np];k++)
				{
					for (s_uint m=0;m<8;m++)
					{
						double x = PositionG[Ele_Np[k]*24+3*m+0];
						double y = PositionG[Ele_Np[k]*24+3*m+1];
						double z = PositionG[Ele_Np[k]*24+3*m+2];
						p <<1,
							x,
							y,
							z,
							x*y,
							x*z,
							y*z,
		                    x*x,
							y*y,
							z*z;
						
						A = A + p * p.transpose();
						b = b + p * stressG[Ele_Np[k]*48+6*m+j];
					}					
				}
				a = A.inverse()*b;
				for (s_uint k=0;k<N_SPR;k++)
				{
					//已知待求节点编号是NP_SPR[k],求它所在的单元和单元中的节点号
					uint N_ele; 
					s_uint Np_ele;
					for(uint Ele = 0; Ele < NUME; Ele++ )
					{   
						for (s_uint N=0;N<8;N++)
						{
							if (Ele_NodeNumber[Ele*8+N]==NP_SPR[k])
							{
								N_ele = Ele;
								Np_ele = N;
								break;
							}
						}
					}
			       double x = PrePositions[N_ele*24+Np_ele*3+0];
				   double y = PrePositions[N_ele*24+Np_ele*3+1];
				   double z = PrePositions[N_ele*24+Np_ele*3+2];				
						p <<1,
							x,
							y,
							z,
							x*y,
							x*z,
							y*z,
		                    x*x,
							y*y,
							z*z;
						MatrixXd tmp;
						tmp =  p.transpose()*a;
						stress_tmp[(NP_SPR[k]-1)*6+j] = stress_tmp[(NP_SPR[k]-1)*6+j] + tmp(0,0);
				}
			}

		}
	}

	for(uint Ele = 0; Ele < NUME; Ele++ )
	{
		for (s_uint N=0;N<8;N++)
		{
			for (s_uint j=0;j<6;j++)
			{
				stress_SPR[Ele*48+N*6+j] = stress_tmp[(Ele_NodeNumber[Ele*8+N]-1)*6+j]/CNT[Ele_NodeNumber[Ele*8+N]-1];
			}
		}
	}
}
