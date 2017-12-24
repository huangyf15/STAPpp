SUBROUTINE shell
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To set up storage and call the truss element subroutine         .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: NUME, NUMMAT, MM, N101, N102, N103, N104, N105, N106,N107,N108,N109

  NUME = NPAR(2)
  NUMMAT = NPAR(3)

! Allocate storage for element group data
  IF (IND == 1) THEN
      MM = 3*NUMMAT*ITWO + 25*NUME + 8*NUME*ITWO
      CALL MEMALLOC(11,"ELEGP",MM,1)
  END IF

  NFIRST=NP(11)   ! Pointer to the first entry in the element group data array
                  ! in the unit of single precision (corresponding to A)

! Calculate the pointer to the arrays in the element group data
! N101: E(NUMMAT)
! N102: miu(NUMMAT)
! N103: h(NUMMAT)
! N104: RHO(NUMMAT)
! N105: LM(20,NUME)
! N106: XY(8,NUME)
! N107: MATP(NUME)
! N108: ELID(4,NUME)
  N101=NFIRST
  N102=N101+NUMMAT*ITWO
  N103=N102+NUMMAT*ITWO  
  N104=N103+NUMMAT*ITWO
  N105=N104+NUMMAT*ITWO
  N106=N105+20*NUME
  N107=N106+8*NUME*ITWO
  N108=N107+NUME
  N109=N108+NUME*4
  NLAST=N109

  MIDEST=NLAST - NFIRST

  CALL shellRUSS (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)),DA(NP(4)),IA(NP(5)),   &
       A(N101),A(N102),A(N103),A(N104),A(N105),A(N106),A(N107),A(N108))

  RETURN

END SUBROUTINE shell


SUBROUTINE shellRUSS (ID,X,Y,Z,U,MHT,E,miu,h,rou,LM,XY,MATP,ELID)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   TRUSS element subroutine                                        .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: ID(6,NUMNP),LM(20,NPAR(2)),MATP(NPAR(2)),MHT(NEQ),ELID(4,NPAR(2)),NNUM(NUMNP)
  REAL(8) :: X(NUMNP),Y(NUMNP),Z(NUMNP),E(NPAR(3)),miu(NPAR(3)),h(NPAR(3)),  &
             XY(8,NPAR(2)),U(NEQ),rou(NPAR(3)),FORCE(3,NUMNP)
  REAL(8) :: S(20,20)
  REAL(8) ::R(NEQ)

  INTEGER :: NPAR1, NUME, NUMMAT, ND, I1, I2, I3, I4, L, N,I,II
  INTEGER :: MTYPE, IPRINT
  REAL(8) :: XL2, XL, SQRT, XX, YY, STR, P
  

  
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  
  INTEGER :: J1,J2,ZZS1,J3
  REAL(8) :: kesai(4),yita(4),kesai0(2),yita0(2),D0(3,3),A1,B1,C1,D1,de(20,1),sigma1(3,1),sigma2(3,1),HH,aa,bb,Be(3,20),Be1(3,20)
  kesai(1)=-1.D0
  kesai(2)=1.D0
  kesai(3)=1.D0
  kesai(4)=-1.D0
  yita(1)=-1.D0
  yita(2)=-1.D0
  yita(3)=1.D0
  yita(4)=1.D0
  kesai0(1)=-0.5773502692
  kesai0(2)=0.5773502692
  yita0(1)=-0.5773502692
  yita0(2)=0.5773502692
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  

  NPAR1  = NPAR(1)
  NUME   = NPAR(2)
  NUMMAT = NPAR(3) 

  ND=20

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
        READ (IIN,'(I10,5F20.0)') N,E(N),miu(N),h(N),rou(N)  ! Read material information
        WRITE (IOUT,"(I5,4X,E12.5,2X,E14.6,2X,E14.6,2X,E14.6,2X,E14.6,2X,E14.6)") N,E(N),miu(N),h(N),rou(N)  
     END DO

     WRITE (IOUT,"(//,' E L E M E N T   I N F O R M A T I O N',//,  &
                      ' ELEMENT     NODE     NODE     NODE     NODE       MATERIAL',/,   &
                      ' NUMBER-N      I        J        J        J       SET NUMBER')")

     N=0
     DO WHILE (N .NE. NUME)
        READ (IIN,'(7I10)') N,I1,I2,I3,I4,MTYPE  ! Read in element information

        ELID(:,N) = (/I1,I2,I3,I4/)
        
!       Save element information
        XY(1,N)=X(I1)  
        XY(2,N)=Y(I1)

        XY(3,N)=X(I2)  
        XY(4,N)=Y(I2)
        
        XY(5,N)=X(I3)  
        XY(6,N)=Y(I3)       
        
        XY(7,N)=X(I4)  
        XY(8,N)=Y(I4)


        MATP(N)=MTYPE  ! Material type
        DO L=1,20
           LM(L,N)=0
        END DO

        DO L=1,5
           LM(L,N)=ID(L,I1)     ! Connectivity matrix
           LM(L+5,N)=ID(L,I2)
           LM(L+10,N)=ID(L,I3)
           LM(L+15,N)=ID(L,I4)
        END DO

!       Update column heights and bandwidth
        CALL COLHT (MHT,ND,LM(1,N))   

        WRITE (IOUT,"(I5,6X,I5,4X,I5,4X,I5,4X,I5,7X,I5)") N,I1,I2,I3,I4,MTYPE

     END DO

     CALL POSTPROCESS(3,4,3,ID,X,Y,Z,U,ELID,FORCE)
     
     RETURN

! Assemble stucture stiffness matrix
  ELSE IF (IND .EQ. 2) THEN
      
      REWIND ILOAD
      READ (ILOAD) R
      

     DO N=1,NUME
        MTYPE=MATP(N)
        
        aa=(abs(XY(1,N)-XY(3,N))+abs(XY(3,N)-XY(5,N))+abs(XY(5,N)-XY(7,N))+abs(XY(7,N)-XY(1,N)))/8.D0
        bb=(abs(XY(2,N)-XY(4,N))+abs(XY(4,N)-XY(6,N))+abs(XY(6,N)-XY(8,N))+abs(XY(8,N)-XY(2,N)))/8.D0
        
        S=0;
        HH=E(MTYPE)*h(MTYPE)**3/(12.D0*(1.D0-miu(MTYPE)*miu(MTYPE)))/(60.D0*aa*bb)
        
        D0(1,1)=E(MTYPE)/(1.D0-miu(MTYPE)*miu(MTYPE))
        D0(2,2)=D0(1,1)
        D0(1,2)=E(MTYPE)*miu(MTYPE)/(1.D0-miu(MTYPE)*miu(MTYPE))
        D0(2,1)=D0(1,2)
        D0(3,3)=E(MTYPE)/(2.D0+2.D0*miu(MTYPE))
        
!        DO J1=1,4
!           DO J2=1,4
!               S(5*J1-2,5*J2-2)=3*HH*(15*((bb*bb)*kesai(J1)*kesai(J2)/(aa*aa)+(aa*aa)*yita(J1)*yita(J2)/(bb*bb))+(14-4*miu(MTYPE)+5*(bb*bb)/(aa*aa)+5*(aa*aa)/(bb*bb))*kesai(J1)*kesai(J2)*yita(J1)*yita(J2))
!               S(5*J1-2,5*J2-1)=-3*HH*bb*((2+3*miu(MTYPE)+5*(aa*aa)/(bb*bb))*kesai(J1)*kesai(J2)*yita(J1)+15*(aa*aa)/(bb*bb)*yita(J1)+5*miu(MTYPE)*kesai(J1)*kesai(J2)*yita(J2))
!               S(5*J1-2,5*J2)=3*HH*aa*((2+3*miu(MTYPE)+5*(bb*bb)/(aa*aa))*kesai(J1)*yita(J1)*yita(J2)+15*(bb*bb)/(aa*aa)*kesai(J1)+5*miu(MTYPE)*kesai(J2)*yita(J1)*yita(J2))
!               S(5*J1-1,5*J2-2)=-3*HH*bb*((2+3*miu(MTYPE)+5*(aa*aa)/(bb*bb))*kesai(J1)*kesai(J2)*yita(J2)+15*(aa*aa)/(bb*bb)*yita(J2)+5*miu(MTYPE)*kesai(J1)*kesai(J2)*yita(J1))
!               S(5*J1-1,5*J2-1)=HH*bb*bb*(2*(1-miu(MTYPE))*kesai(J1)*kesai(J2)*(3+5*yita(J1)*yita(J2))+5*(aa*aa)/(bb*bb)*(3+kesai(J1)*kesai(J2))*(3+yita(J1)*yita(J2)))
!               S(5*J1-1,5*J2)=-15*HH*miu(MTYPE)*aa*bb*(kesai(J1)+kesai(J2))*(yita(J1)+yita(J2))
!               S(5*J1,5*J2-2)=3*HH*aa*((2+3*miu(MTYPE)+5*(bb*bb)/(aa*aa))*kesai(J2)*yita(J1)*yita(J2)+15*(bb*bb)/(aa*aa)*kesai(J2)+5*miu(MTYPE)*kesai(J1)*yita(J1)*yita(J2))
!               S(5*J1,5*J2-1)=-15*HH*miu(MTYPE)*aa*bb*(kesai(J1)+kesai(J2))*(yita(J1)+yita(J2))
!               S(5*J1,5*J2)=HH*aa*aa*(2*(1-miu(MTYPE))*yita(J1)*yita(J2)*(3+5**kesai(J1)*kesai(J2))+5*(bb*bb)/(aa*aa)*(3+kesai(J1)*kesai(J2))*(3+yita(J1)*yita(J2)))
!            enddo
!        enddo
   
        DO J1=1,2
            DO J2=1,2
                Be=0.D0
                Be(1,1)=(yita0(J1)-1.D0)/(4.D0*aa)
                Be(1,6)=-(yita0(J1)-1.D0)/(4.D0*aa)
                Be(1,11)=(yita0(J1)+1.D0)/(4.D0*aa)
                Be(1,16)=-(yita0(J1)+1.D0)/(4.D0*aa)
                Be(2,2)=(kesai0(J2)-1.D0)/(4*bb)
                Be(2,7)=-(kesai0(J2)+1.D0)/(4*bb)
                Be(2,12)=(kesai0(J2)+1.D0)/(4*bb)
                Be(2,17)=(-kesai0(J2)+1.D0)/(4*bb)
                Be(3,1)=Be(2,2)
                Be(3,2)=Be(1,1)
                Be(3,6)=Be(2,7)
                Be(3,7)=Be(1,6)
                Be(3,11)=Be(2,12)
                Be(3,12)=Be(1,11)
                Be(3,16)=Be(2,17)
                Be(3,17)=Be(1,16)
                
                Be1=0.D0
                DO J3=1,4
                    Be1(1,5*J3-2)=(3.D0/(4.D0*aa*aa))*kesai(J3)*kesai0(J2)*(1.D0+yita(J3)*yita0(J1))
                    Be1(1,5*J3)=(1.D0/(4.D0*aa))*(1+yita(J3)*yita0(J1))*(3*kesai0(J2)+kesai(J3))
                    Be1(2,5*J3-2)=(3.D0/(4.D0*bb*bb))*yita(J3)*yita0(J1)*(1.D0+kesai(J3)*kesai0(J2))
                    Be1(2,5*J3-1)=-(1.D0/(4.D0*bb))*(1.D0+kesai(J3)*kesai0(J2))*(3.D0*yita0(J1)+yita(J3))
                    Be1(3,5*J3-2)=(1.D0/(4.D0*aa*bb))*kesai(J3)*yita(J3)*(3.D0*kesai0(J2)*kesai0(J2)+3.D0*yita0(J1)*yita0(J1)-4.D0)
                    Be1(3,5*J3-1)=(1.D0/(4.D0*aa))*kesai(J3)*(1.D0-2.D0*yita(J3)*yita0(J1)-3*yita0(J1)*yita0(J1))
                    Be1(3,5*J3)=-(1.D0/(4.D0*bb))*yita(J3)*(1.D0-2.D0*kesai(J3)*kesai0(J2)-3*kesai0(J2)*kesai0(J2)) 
                enddo
                
                S=S+matmul(matmul(transpose(Be),D0*h(MTYPE)),Be)*(aa*bb)+matmul(matmul(transpose(Be1),D0),Be1)*(aa*bb)*h(MTYPE)**3/12.D0
            enddo
        enddo
        
     
        CALL ADDBAN (DA(NP(3)),IA(NP(2)),S,LM(1,N),ND)
      
        DO I=1,4
            II=LM(5*I-4,N)
            IF(II>0) R(II)=R(II)+gravity(1)*h(MTYPE)*rou(MTYPE)*aa*bb/4;	
            II=LM(5*I-3,N)	
            IF(II>0) R(II)=R(II)+gravity(2)*h(MTYPE)*rou(MTYPE)*aa*bb/4;	
            II=LM(5*I-2,N)
            IF(II>0) R(II)=R(II)+gravity(3)*h(MTYPE)*rou(MTYPE)*aa*bb/4;	
        ENDDO            

     END DO
     
      REWIND ILOAD
      WRITE (ILOAD) R

     RETURN

! Stress calculations
  ELSE IF (IND .EQ. 3) THEN

      NNUM = 0
      FORCE = 0.D0

     IPRINT=0
     DO N=1,NUME
        IPRINT=IPRINT + 1
        IF (IPRINT.GT.50) IPRINT=1
        IF (IPRINT.EQ.1) WRITE (IOUT,"(//,' S T R E S S  C A L C U L A T I O N S  F O R  ',  &
                                           'E L E M E N T  G R O U P',I4,//,   &
                                           '  ELEMENT',13X,'x',9X,'y',18X,'STRESS',1X,'xx',9X,'yy',9X,'xy',/,'  NUMBER')") NG
        MTYPE=MATP(N)

        aa=(abs(XY(1,N)-XY(3,N))+abs(XY(3,N)-XY(5,N))+abs(XY(5,N)-XY(7,N))+abs(XY(7,N)-XY(1,N)))/8.D0
        bb=(abs(XY(2,N)-XY(4,N))+abs(XY(4,N)-XY(6,N))+abs(XY(6,N)-XY(8,N))+abs(XY(8,N)-XY(2,N)))/8.D0
        
        S=0;
        HH=E(MTYPE)*h(MTYPE)**3/(12.D0*(1.D0-miu(MTYPE)*miu(MTYPE)))/(60.D0*aa*bb)
        
        D0(1,1)=E(MTYPE)/(1.D0-miu(MTYPE)*miu(MTYPE))
        D0(2,2)=D0(1,1)
        D0(1,2)=E(MTYPE)*miu(MTYPE)/(1.D0-miu(MTYPE)*miu(MTYPE))
        D0(2,1)=D0(1,2)
        D0(3,3)=E(MTYPE)/(2.D0+2.D0*miu(MTYPE))
        
        DO L=1,20
           I1=LM(L,N)
             IF (I1 == 0) THEN
             de(L,1)=0.D0
             ELSE 
             de(L,1)=U(I1)
             END IF
           
        END DO

        
        DO J1=1,2
            DO J2=1,2
                Be=0.D0
                Be(1,1)=(yita0(J1)-1.D0)/(4.D0*aa)
                Be(1,6)=-(yita0(J1)-1.D0)/(4.D0*aa)
                Be(1,11)=(yita0(J1)+1.D0)/(4.D0*aa)
                Be(1,16)=-(yita0(J1)+1.D0)/(4.D0*aa)
                Be(2,2)=(kesai0(J2)-1.D0)/(4.D0*bb)
                Be(2,7)=-(kesai0(J2)+1.D0)/(4.D0*bb)
                Be(2,12)=(kesai0(J2)+1.D0)/(4.D0*bb)
                Be(2,17)=(-kesai0(J2)+1.D0)/(4.D0*bb)
                Be(3,1)=Be(2,2)
                Be(3,2)=Be(1,1)
                Be(3,6)=Be(2,7)
                Be(3,7)=Be(1,6)
                Be(3,11)=Be(2,12)
                Be(3,12)=Be(1,11)
                Be(3,16)=Be(2,17)
                Be(3,17)=Be(1,16)
                
                Be1=0.D0
                DO J3=1,4
                    Be1(1,5*J3-2)=(3.D0/(4.D0*aa*aa))*kesai(J3)*kesai0(J2)*(1.D0+yita(J3)*yita0(J1))
                    Be1(1,5*J3)=(1.D0/(4.D0*aa))*(1.D0+yita(J3)*yita0(J1))*(3.D0*kesai0(J2)+kesai(J3))
                    Be1(2,5*J3-2)=(3.D0/(4.D0*bb*bb))*yita(J3)*yita0(J1)*(1.D0+kesai(J3)*kesai0(J2))
                    Be1(2,5*J3-1)=-(1.D0/(4.D0*bb))*(1.D0+kesai(J3)*kesai0(J2))*(3.D0*yita0(J1)+yita(J3))
                    Be1(3,5*J3-2)=(1.D0/(4.D0*aa*bb))*kesai(J3)*yita(J3)*(3.D0*kesai0(J2)*kesai0(J2)+3.D0*yita0(J1)*yita0(J1)-4.D0)
                    Be1(3,5*J3-1)=(1.D0/(4.D0*aa))*kesai(J3)*(1.D0-2.D0*yita(J3)*yita0(J1)-3.D0*yita0(J1)*yita0(J1))
                    Be1(3,5*J3)=-(1.D0/(4.D0*bb))*yita(J3)*(1.D0-2.D0*kesai(J3)*kesai0(J2)-3.D0*kesai0(J2)*kesai0(J2)) 
                enddo
                sigma1=matmul(-h(MTYPE)*D0/2,matmul(Be1,de))+matmul(D0*h(MTYPE),matmul(Be,de))
                sigma2=matmul(D0*h(MTYPE),matmul(Be,de))
                WRITE (IOUT,"(1X,I5,11X,E9.2,4X,E9.2,4X,E9.2,4X,E9.2,4X,E9.2,4X,E9.2)") N,sigma1(1,1),sigma1(2,1),sigma1(3,1),sigma2(1,1),sigma2(2,1),sigma2(3,1)
        
        DO I=1,3
            FORCE(:,ELID(I,N)) = FORCE(:,ELID(I,N)) + SIGMA2(:,1)
        END DO
        NNUM(ELID(:,N)) = NNUM(ELID(:,N)) + 1
           
            enddo
        enddo
        
   

  !      WRITE (IOUT,"(1X,I5,11X,E13.6,4X,E13.6)") N,P,STR
     END DO
     
    ! For postprocess
     DO I = 1,NUMNP
         IF(NNUM(I) .NE. 0) FORCE(:,I) = FORCE(:,I)/NNUM(I)
     END DO
    CALL POSTPROCESS(3,4,3,ID,X,Y,(/(0.D0, I=1,NUMNP)/),U,ELID,FORCE)

  ELSE 
     STOP " ERROR"
  END IF

END SUBROUTINE shellRUSS