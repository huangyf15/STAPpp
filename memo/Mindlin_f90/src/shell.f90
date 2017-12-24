SUBROUTINE shell
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To set up storage and call the truss element subroutine         .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: NUME, NUMMAT, MM, N101, N102, N103, N104, N105, N106,N107

  NUME = NPAR(2)
  NUMMAT = NPAR(3)

! Allocate storage for element group data
  IF (IND == 1) THEN
      MM = 3*NUMMAT*ITWO + 41*NUME + 16*NUME*ITWO
      CALL MEMALLOC(11,"ELEGP",MM,1)
  END IF

  NFIRST=NP(11)   ! Pointer to the first entry in the element group data array
                  ! in the unit of single precision (corresponding to A)

! Calculate the pointer to the arrays in the element group data
! N101: E(NUMMAT)
! N102: miu(NUMMAT)
! N103: h(NUMMAT)
! N104: LM(40,NUME)
! N105: XY(16,NUME)
! N106: MTAP(NUME)
  N101=NFIRST
  N102=N101+NUMMAT*ITWO
  N103=N102+NUMMAT*ITWO  
  N104=N103+NUMMAT*ITWO
  N105=N104+40*NUME
  N106=N105+16*NUME*ITWO
  N107=N106+NUME
  NLAST=N107

  MIDEST=NLAST - NFIRST

  CALL shellRUSS (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)),IA(NP(5)),   &
       A(N101),A(N102),A(N103),A(N104),A(N105),A(N106))

  RETURN

END SUBROUTINE shell


SUBROUTINE shellRUSS (ID,X,Y,U,MHT,E,miu,h,LM,XY,MATP)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   TRUSS element subroutine                                        .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: ID(6,NUMNP),LM(40,NPAR(2)),MATP(NPAR(2)),MHT(NEQ)
  REAL(8) :: X(NUMNP),Y(NUMNP),E(NPAR(3)),miu(NPAR(3)),h(NPAR(3)),  &
             XY(16,NPAR(2)),U(NEQ)
  REAL(8) :: S(40,40)

  INTEGER :: NPAR1, NUME, NUMMAT, ND, I1, I2, I3, I4,I5,I6,I7,I8, L, N,I
  INTEGER :: MTYPE, IPRINT
  REAL(8) :: XL2, XL, SQRT, XX, YY, STR, P
  

  
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  
  INTEGER :: J1,J2,ZZS1
  REAL(8) :: kesai(3),yita(3),Be1(3,40),D0(3,3),Be2(2,40),Be3(3,40),A1,B1,C1,D1,de(40,1),sigma1(3,1),sigma2(3,1),sigma3(2,1),Nx(2,8),NN(2,8),JJ(2,2),WWW(3),NN0(1,8)
  kesai(1)=-SQRT(0.6D0)
  kesai(2)=0
  kesai(3)=SQRT(0.6D0)
  yita(1)=-SQRT(0.6D0)
  yita(2)=0.D0
  yita(3)=SQRT(0.6D0)
  WWW(1)=(5.D0/9.D0)
  WWW(2)=(8.D0/9.D0)
  WWW(3)=(5.D0/9.D0)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  

  NPAR1  = NPAR(1)
  NUME   = NPAR(2)
  NUMMAT = NPAR(3) 

  ND=40

! Read and generate element information
  IF (IND .EQ. 1) THEN

     WRITE (IOUT,"(' E L E M E N T   D E F I N I T I O N',//,  &
                   ' ELEMENT TYPE ',13(' .'),'( NPAR(1) ) . . =',I5,/,   &
                   '     EQ.1, 4Q  ELEMENTS',/,      &
                   '     EQ.2, ELEMENTS CURRENTLY',/,  &
                   '     EQ.3, NOT AVAILABLE',//,      &
                   ' NUMBER OF ELEMENTS.',10(' .'),'( NPAR(2) ) . . =',I5,/)") NPAR1,NUME

     IF (NUMMAT.EQ.0) NUMMAT=1

     WRITE (IOUT,"(' M A T E R I A L   D E F I N I T I O N',//,  &
                   ' NUMBER OF DIFFERENT SETS OF MATERIAL',/,  &
                   ' AND CROSS-SECTIONAL  CONSTANTS ',         &
                   4 (' .'),'( NPAR(3) ) . . =',I5,/)") NUMMAT

     WRITE (IOUT,"('  SET       YOUNG''S     CROSS-SECTIONAL',/,  &
                   ' NUMBER     MODULUS',10X,'miu ',/,  &
                   15 X,'E',14X,'A')")

     DO I=1,NUMMAT
        READ (IIN,'(I5,3F10.0)') N,E(N),miu(N),h(N)  ! Read material information
        WRITE (IOUT,"(I5,4X,E12.5,2X,E14.6,2X,E14.6)") N,E(N),miu(N),h(N)
     END DO

     WRITE (IOUT,"(//,' E L E M E N T   I N F O R M A T I O N',//,  &
                      ' ELEMENT     NODE     NODE     NODE     NODE       MATERIAL',/,   &
                      ' NUMBER-N      I        J        J        J       SET NUMBER')")

     N=0
     DO WHILE (N .NE. NUME)
        READ (IIN,'(11I5)') N,I1,I2,I3,I4,I5,I6,I7,I8,MTYPE  ! Read in element information

!       Save element information
        XY(1,N)=X(I1)  
        XY(2,N)=Y(I1)

        XY(3,N)=X(I2)  
        XY(4,N)=Y(I2)
        
        XY(5,N)=X(I3)  
        XY(6,N)=Y(I3)       
        
        XY(7,N)=X(I4)  
        XY(8,N)=Y(I4)
        
        XY(9,N)=X(I5)  
        XY(10,N)=Y(I5)

        XY(11,N)=X(I6)  
        XY(12,N)=Y(I6)
        
        XY(13,N)=X(I7)  
        XY(14,N)=Y(I7)       
        
        XY(15,N)=X(I8)  
        XY(16,N)=Y(I8)


        MATP(N)=MTYPE  ! Material type
        DO L=1,40
           LM(L,N)=0
        END DO

        DO L=1,5
           LM(L,N)=ID(L,I1)     ! Connectivity matrix
           LM(L+5,N)=ID(L,I2)
           LM(L+10,N)=ID(L,I3)
           LM(L+15,N)=ID(L,I4)
           LM(L+20,N)=ID(L,I5)
           LM(L+25,N)=ID(L,I6)
           LM(L+30,N)=ID(L,I7)
           LM(L+35,N)=ID(L,I8)
        END DO

!       Update column heights and bandwidth
        CALL COLHT (MHT,ND,LM(1,N))   

        WRITE (IOUT,"(I5,6X,I5,4X,I5,4X,I5,4X,I5,4X,I5,4X,I5,4X,I5,4X,I5,7X,I5)") N,I1,I2,I3,I4,I5,I6,I7,I8,MTYPE

     END DO

     RETURN

! Assemble stucture stiffness matrix
  ELSE IF (IND .EQ. 2) THEN

     DO N=1,NUME
        MTYPE=MATP(N)
        
        S=0;
        D0(1,1)=E(MTYPE)/(1-miu(MTYPE)*miu(MTYPE))
        D0(2,2)=D0(1,1)
        D0(1,2)=E(MTYPE)*miu(MTYPE)/(1-miu(MTYPE)*miu(MTYPE))
        D0(2,1)=D0(1,2)
        D0(3,3)=E(MTYPE)/(2+2*miu(MTYPE))
        
        

        DO J1=1,3
            DO J2=1,3
                
            NN0(1,1)=(1-kesai(J2))*(1-yita(J1))*(-kesai(J2)-yita(J1)-1)/4
            NN0(1,2)=(1+kesai(J2))*(1-yita(J1))*(kesai(J2)-yita(J1)-1)/4
            NN0(1,3)=(1+kesai(J2))*(1+yita(J1))*(kesai(J2)+yita(J1)-1)/4
            NN0(1,4)=(1-kesai(J2))*(1+yita(J1))*(-kesai(J2)+yita(J1)-1)/4
            NN0(1,5)=(1-kesai(J2)*kesai(J2))*(1-yita(J1))/2
            NN0(1,6)=(1-yita(J1)*yita(J1))*(1+kesai(J2))/2
            NN0(1,7)=(1-kesai(J2)*kesai(J2))*(1+yita(J1))/2
            NN0(1,8)=(1-yita(J1)*yita(J1))*(1-kesai(J2))/2
            
              NN(1,1)=(2.D0*kesai(J2)+yita(J1))*(1.D0-yita(J1))/4.D0
              NN(2,1)=(2.D0*yita(J1)+kesai(J2))*(1.D0-kesai(J2))/4.D0
              NN(1,2)=(2.D0*kesai(J2)-yita(J1))*(1.D0-yita(J1))/4.D0
              NN(2,2)=(2.D0*yita(J1)-kesai(J2))*(1.D0+kesai(J2))/4.D0
              NN(1,3)=(2.D0*kesai(J2)+yita(J1))*(1.D0+yita(J1))/4.D0
              NN(2,3)=(2.D0*yita(J1)+kesai(J2))*(1.D0+kesai(J2))/4.D0
              NN(1,4)=(2.D0*kesai(J2)-yita(J1))*(1.D0+yita(J1))/4.D0
              NN(2,4)=(2.D0*yita(J1)-kesai(J2))*(1.D0-kesai(J2))/4.D0
              NN(1,5)=-kesai(J2)*(1.D0-yita(J1))
              NN(2,5)=-(1.D0-kesai(J2)*kesai(J2))/2.D0
              NN(1,6)=(1.D0-yita(J1)*yita(J1))/2.D0
              NN(2,6)=-yita(J1)*(1.D0+kesai(J2))
              NN(1,7)=-kesai(J2)*(1.D0+yita(J1))
              NN(2,7)=(1.D0-kesai(J2)*kesai(J2))/2.D0
              NN(1,8)=-(1.D0-yita(J1)*yita(J1))/2.D0
              NN(2,8)=-yita(J1)*(1.D0-kesai(J2))
                
                A1=NN(1,1)*XY(1,N)+NN(1,2)*XY(3,N)+NN(1,3)*XY(5,N)+NN(1,4)*XY(7,N)+NN(1,5)*XY(9,N)+NN(1,6)*XY(11,N)+NN(1,7)*XY(13,N)+NN(1,8)*XY(15,N)
                B1=NN(1,1)*XY(2,N)+NN(1,2)*XY(4,N)+NN(1,3)*XY(6,N)+NN(1,4)*XY(8,N)+NN(1,5)*XY(10,N)+NN(1,6)*XY(12,N)+NN(1,7)*XY(14,N)+NN(1,8)*XY(16,N)
                C1=NN(2,1)*XY(1,N)+NN(2,2)*XY(3,N)+NN(2,3)*XY(5,N)+NN(2,4)*XY(7,N)+NN(2,5)*XY(9,N)+NN(2,6)*XY(11,N)+NN(2,7)*XY(13,N)+NN(2,8)*XY(15,N)
                D1=NN(2,1)*XY(2,N)+NN(2,2)*XY(4,N)+NN(2,3)*XY(6,N)+NN(2,4)*XY(8,N)+NN(2,5)*XY(10,N)+NN(2,6)*XY(12,N)+NN(2,7)*XY(14,N)+NN(2,8)*XY(16,N)
                
                JJ(1,1)=D1/(A1*D1-B1*C1)
                JJ(1,2)=-B1/(A1*D1-B1*C1)
                JJ(2,1)=-C1/(A1*D1-B1*C1)
                JJ(2,2)=A1/(A1*D1-B1*C1)
                
                Nx=matmul(JJ,NN)
                
                
                Be1=0
                Be3=0
                Be2=0
                
                DO ZZS1=1,8
                
                Be1(1,5*ZZS1-1)=Nx(1,ZZS1)               
                Be1(2,5*ZZS1)=Nx(2,ZZS1)
                Be1(3,5*ZZS1-1)=Nx(2,ZZS1)
                Be1(3,5*ZZS1)=Nx(1,ZZS1)
                Be3(1,5*ZZS1-4)=Nx(1,ZZS1)              
                Be3(2,5*ZZS1-3)=Nx(2,ZZS1)
                Be3(3,5*ZZS1-3)=Nx(1,ZZS1)
                Be3(3,5*ZZS1-4)=Nx(2,ZZS1)
                Be2(1,5*ZZS1-2)=Nx(1,ZZS1)
                Be2(2,5*ZZS1-2)=Nx(2,ZZS1)
                Be2(1,5*ZZS1-1)=-NN0(1,ZZS1)
                Be2(2,5*ZZS1)=-NN0(1,ZZS1)
                
                ENDDO

                
                S=S+WWW(J1)*WWW(J2)*(matmul(matmul(transpose(Be1),D0*h(MTYPE)**3/12),Be1)*abs(A1*D1-B1*C1)+matmul(matmul(transpose(Be3),D0*h(MTYPE)),Be3)*abs(A1*D1-B1*C1)+matmul(transpose(Be2),E(MTYPE)/(2+2*miu(MTYPE))*h(MTYPE)*Be2)*abs(A1*D1-B1*C1)*1.2)

            enddo
        enddo
        
        CALL ADDBAN (DA(NP(3)),IA(NP(2)),S,LM(1,N),ND)

     END DO

     RETURN

! Stress calculations
  ELSE IF (IND .EQ. 3) THEN

     IPRINT=0
     DO N=1,NUME
        IPRINT=IPRINT + 1
        IF (IPRINT.GT.50) IPRINT=1
        IF (IPRINT.EQ.1) WRITE (IOUT,"(//,' S T R E S S  C A L C U L A T I O N S  F O R  ',  &
                                           'E L E M E N T  G R O U P',I4,//,   &
                                           '  ELEMENT',13X,'x',9X,'y',18X,'STRESS',1X,'xx',9X,'yy',9X,'xy',/,'  NUMBER')") NG
        MTYPE=MATP(N)

        D0(1,1)=E(MTYPE)/(1-miu(MTYPE)*miu(MTYPE))
        D0(2,2)=D0(1,1)
        D0(1,2)=E(MTYPE)*miu(MTYPE)/(1-miu(MTYPE)*miu(MTYPE))
        D0(2,1)=D0(1,2)
        D0(3,3)=E(MTYPE)/(2+2*miu(MTYPE))
        
        DO L=1,40
           I1=LM(L,N)
             IF (I1 == 0) THEN
             de(L,1)=0
             ELSE 
             de(L,1)=U(I1)
             END IF
           
        END DO
              

       DO J1=1,3
            DO J2=1,3
                
            NN0(1,1)=(1-kesai(J2))*(1-yita(J1))*(-kesai(J2)-yita(J1)-1)/4
            NN0(1,2)=(1+kesai(J2))*(1-yita(J1))*(kesai(J2)-yita(J1)-1)/4
            NN0(1,3)=(1+kesai(J2))*(1+yita(J1))*(kesai(J2)+yita(J1)-1)/4
            NN0(1,4)=(1-kesai(J2))*(1+yita(J1))*(-kesai(J2)+yita(J1)-1)/4
            NN0(1,5)=(1-kesai(J2)*kesai(J2))*(1-yita(J1))/2
            NN0(1,6)=(1-yita(J1)*yita(J1))*(1+kesai(J2))/2
            NN0(1,7)=(1-kesai(J2)*kesai(J2))*(1+yita(J1))/2
            NN0(1,8)=(1-yita(J1)*yita(J1))*(1-kesai(J2))/2
            
              NN(1,1)=(2.D0*kesai(J2)+yita(J1))*(1.D0-yita(J1))/4.D0
              NN(2,1)=(2.D0*yita(J1)+kesai(J2))*(1.D0-kesai(J2))/4.D0
              NN(1,2)=(2.D0*kesai(J2)-yita(J1))*(1.D0-yita(J1))/4.D0
              NN(2,2)=(2.D0*yita(J1)-kesai(J2))*(1.D0+kesai(J2))/4.D0
              NN(1,3)=(2.D0*kesai(J2)+yita(J1))*(1.D0+yita(J1))/4.D0
              NN(2,3)=(2.D0*yita(J1)+kesai(J2))*(1.D0+kesai(J2))/4.D0
              NN(1,4)=(2.D0*kesai(J2)-yita(J1))*(1.D0+yita(J1))/4.D0
              NN(2,4)=(2.D0*yita(J1)-kesai(J2))*(1.D0-kesai(J2))/4.D0
              NN(1,5)=-kesai(J2)*(1.D0-yita(J1))
              NN(2,5)=-(1.D0-kesai(J2)*kesai(J2))/2.D0
              NN(1,6)=(1.D0-yita(J1)*yita(J1))/2.D0
              NN(2,6)=-yita(J1)*(1.D0+kesai(J2))
              NN(1,7)=-kesai(J2)*(1.D0+yita(J1))
              NN(2,7)=(1.D0-kesai(J2)*kesai(J2))/2.D0
              NN(1,8)=-(1.D0-yita(J1)*yita(J1))/2.D0
              NN(2,8)=-yita(J1)*(1.D0-kesai(J2))
                
                A1=NN(1,1)*XY(1,N)+NN(1,2)*XY(3,N)+NN(1,3)*XY(5,N)+NN(1,4)*XY(7,N)+NN(1,5)*XY(9,N)+NN(1,6)*XY(11,N)+NN(1,7)*XY(13,N)+NN(1,8)*XY(15,N)
                B1=NN(1,1)*XY(2,N)+NN(1,2)*XY(4,N)+NN(1,3)*XY(6,N)+NN(1,4)*XY(8,N)+NN(1,5)*XY(10,N)+NN(1,6)*XY(12,N)+NN(1,7)*XY(14,N)+NN(1,8)*XY(16,N)
                C1=NN(2,1)*XY(1,N)+NN(2,2)*XY(3,N)+NN(2,3)*XY(5,N)+NN(2,4)*XY(7,N)+NN(2,5)*XY(9,N)+NN(2,6)*XY(11,N)+NN(2,7)*XY(13,N)+NN(2,8)*XY(15,N)
                D1=NN(2,1)*XY(2,N)+NN(2,2)*XY(4,N)+NN(2,3)*XY(6,N)+NN(2,4)*XY(8,N)+NN(2,5)*XY(10,N)+NN(2,6)*XY(12,N)+NN(2,7)*XY(14,N)+NN(2,8)*XY(16,N)
                
                JJ(1,1)=D1/(A1*D1-B1*C1)
                JJ(1,2)=-B1/(A1*D1-B1*C1)
                JJ(2,1)=-C1/(A1*D1-B1*C1)
                JJ(2,2)=A1/(A1*D1-B1*C1)
                
                Nx=matmul(JJ,NN)
                
                
                Be1=0
                Be3=0
                Be2=0
                
                DO ZZS1=1,8
                
                Be1(1,5*ZZS1-1)=Nx(1,ZZS1)               
                Be1(2,5*ZZS1)=Nx(2,ZZS1)
                Be1(3,5*ZZS1-1)=Nx(2,ZZS1)
                Be1(3,5*ZZS1)=Nx(1,ZZS1)
                Be3(1,5*ZZS1-4)=Nx(1,ZZS1)              
                Be3(2,5*ZZS1-3)=Nx(2,ZZS1)
                Be3(3,5*ZZS1-3)=Nx(1,ZZS1)
                Be3(3,5*ZZS1-4)=Nx(2,ZZS1)
                Be2(1,5*ZZS1-2)=Nx(1,ZZS1)
                Be2(2,5*ZZS1-2)=Nx(2,ZZS1)
                Be2(1,5*ZZS1-1)=-NN0(1,ZZS1)
                Be2(2,5*ZZS1)=-NN0(1,ZZS1)
                
                ENDDO

                
               
                sigma1=matmul(-h(MTYPE)*D0/2,matmul(Be1,de))+matmul(D0,matmul(Be3,de))
                sigma2=matmul(D0,matmul(Be3,de))
                sigma3=matmul(E(MTYPE)/(2+2*miu(MTYPE))*Be2,de)
                
                WRITE (IOUT,"(1X,I5,11X,E9.2,4X,E9.2,4X,E9.2,4X,E9.2,4X,E9.2,4X,E9.2,4X,E9.2,4X,E9.2)") N,sigma1(1,1),sigma1(2,1),sigma1(3,1),sigma2(1,1),sigma2(2,1),sigma2(3,1),sigma3(1,1),sigma3(2,1)
            enddo
        enddo
                




  !      WRITE (IOUT,"(1X,I5,11X,E13.6,4X,E13.6)") N,P,STR
     END DO

  ELSE 
     STOP "ERROR"
  END IF

END SUBROUTINE shellRUSS