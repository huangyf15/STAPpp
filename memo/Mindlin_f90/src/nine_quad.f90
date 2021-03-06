SUBROUTINE NINE_QUAD
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                               .
! .   To set up storage and call the 9Q element subroutine        .
! .                                                               .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: NUME, NUMMAT, MM, N101, N102, N103, N104, N105, N106

  NUME = NPAR(2)
  NUMMAT = NPAR(3)

! Allocate storage for element group data
  IF (IND == 1) THEN
      MM = 2*NUMMAT*ITWO + 19*NUME + 18*NUME*ITWO
      CALL MEMALLOC(11,"ELEGP",MM,1)
  END IF

  NFIRST=NP(11)   ! Pointer to the first entry in the element group data array
                  ! in the unit of single precision (corresponding to A)

! Calculate the pointer to the arrays in the element group data
! N101: E(NUMMAT)
! N102: AREA(NUMMAT)
! N103: LM(18,NUME)
! N104: XYZ(18,NUME)
! N105: MTAP(NUME)
  N101=NFIRST
  N102=N101+NUMMAT*ITWO
  N103=N102+NUMMAT*ITWO
  N104=N103+18*NUME
  N105=N104+18*NUME*ITWO
  N106=N105+NUME
  NLAST=N106

  MIDEST=NLAST - NFIRST

  CALL INE_QUAD (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)),IA(NP(5)),   &
       A(N101),A(N102),A(N103),A(N104),A(N105))

  RETURN

<<<<<<< HEAD
END SUBROUTINE NINE_QUAD

=======
END SUBROUTINE NINE_QUAD

>>>>>>> master
SUBROUTINE INE_QUAD(ID,X,Y,U,MHT,E,POISSON,LM,XY,MATP)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   9Q element subroutine                                           .
! .                                                                   .
<<<<<<< HEAD
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: ID(6,NUMNP),LM(18,NPAR(2)),MATP(NPAR(2)),MHT(NEQ)
  REAL(8) :: X(NUMNP),Y(NUMNP),E(NPAR(3)),POISSON(NPAR(3)),  &
             XY(18,NPAR(2)),U(NEQ)
  REAL(8) :: JE(2,2),DETJE,JENI(2,2),DERIVAN(2,9),S(18,18),D(3,3)
  REAL(8) :: XI,ETA,TWOWEIGHT
  REAL(8) :: NXI(3),NETA(3),DNXI(3),DNETA(3)
  INTEGER :: POSITION(3,3)
  
  INTEGER :: NPAR1, NUME, NUMMAT, ND, I(9), J1, J2, J3, J4, L, N
  INTEGER :: MTYPE, IPRINT
  REAL(8) :: POISSONE
  REAL(8) :: STR(3), P(3)
  
  INTEGER,PARAMETER :: NGP=3
=======
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: ID(6,NUMNP),LM(18,NPAR(2)),MATP(NPAR(2)),MHT(NEQ)
  REAL(8) :: X(NUMNP),Y(NUMNP),E(NPAR(3)),POISSON(NPAR(3)),  &
             XY(18,NPAR(2)),U(NEQ)
  REAL(8) :: JE(2,2),DETJE,JENI(2,2),DERIVAN(2,9),S(18,18),D(3,3)
  REAL(8) :: XI,ETA,TWOWEIGHT
  REAL(8) :: NXI(3),NETA(3),DNXI(3),DNETA(3)
  INTEGER :: POSITION(3,3)
  
  INTEGER :: NPAR1, NUME, NUMMAT, ND, I(9), J1, J2, J3, J4, L, N
  INTEGER :: MTYPE, IPRINT
  REAL(8) :: POISSONE
  REAL(8) :: STR(3), P(3)
  
  INTEGER,PARAMETER :: NGP=3
>>>>>>> master
  REAL(8) :: GAUSSP(NGP), WEIGHT(NGP)

  NPAR1  = NPAR(1)
  NUME   = NPAR(2)
  NUMMAT = NPAR(3) 

<<<<<<< HEAD
  ND=18
  
  DATA ((POSITION(J1,J2),J2=1,3),J1=1,3)/1, 4, 8, &
                                         2, 3, 6, &
                                         5, 7, 9/
  
=======
  ND=18
  
  DATA ((POSITION(J1,J2),J2=1,3),J1=1,3)/1, 4, 8, &
                                         2, 3, 6, &
                                         5, 7, 9/
  
>>>>>>> master
  ! Read and generate element information
  IF (IND .EQ. 1) THEN

     WRITE (IOUT,"(' E L E M E N T   D E F I N I T I O N',//,  &
                   ' ELEMENT TYPE ',13(' .'),'( NPAR(1) ) . . =',I5,/,   &
                   '     EQ.1, TRUSS ELEMENTS',/,      &
<<<<<<< HEAD
                   '     EQ.2, FOUR QUAD ELEMENTS',/,  &
                   '     EQ.3, ELEMENTS CURRENTLY',/,  &
                   '     EQ.4, NOT AVAILABLE',//,      &
                   ' NUMBER OF ELEMENTS.',10(' .'),'( NPAR(2) ) . . =',I5,/)") NPAR1,NUME
     
     IF (NUMMAT.EQ.0) NUMMAT=1
     
=======
                   '     EQ.2, FOUR QUAD ELEMENTS',/,  &
                   '     EQ.3, ELEMENTS CURRENTLY',/,  &
                   '     EQ.4, NOT AVAILABLE',//,      &
                   ' NUMBER OF ELEMENTS.',10(' .'),'( NPAR(2) ) . . =',I5,/)") NPAR1,NUME
     
     IF (NUMMAT.EQ.0) NUMMAT=1
     
>>>>>>> master
     WRITE (IOUT,"(' M A T E R I A L   D E F I N I T I O N',//,  &
                   ' NUMBER OF DIFFERENT SETS OF MATERIAL',/,  &
                   ' AND CROSS-SECTIONAL  CONSTANTS ',         &
                   4 (' .'),'( NPAR(3) ) . . =',I5,/)") NUMMAT
<<<<<<< HEAD

     WRITE (IOUT,"('  SET       YOUNG''S        POISSON',/,  &
                   ' NUMBER     MODULUS',9X,'RATIO',/,  &
                   15 X,'E',14X,'v')")
     
=======

     WRITE (IOUT,"('  SET       YOUNG''S        POISSON',/,  &
                   ' NUMBER     MODULUS',9X,'RATIO',/,  &
                   15 X,'E',14X,'v')")
     
>>>>>>> master
     DO J1=1,NUMMAT
        READ (IIN,'(I5,2F10.0)') N,E(N),POISSON(N)  ! Read material information
        WRITE (IOUT,"(I5,4X,E12.5,2X,E14.6)") N,E(N),POISSON(N)
     END DO
<<<<<<< HEAD

=======

>>>>>>> master
     WRITE (IOUT,"(//,' E L E M E N T   I N F O R M A T I O N',//,  &
                      ' ELEMENT',9('     NODE'),'       MATERIAL',/,   &
                      ' NUMBER-N','     I',I1,8('       I',I1),'       SET NUMBER')")(J1,J1=1,9)

<<<<<<< HEAD
     N=0
     DO WHILE (N .NE. NUME)
        READ (IIN,'(12I5)') N,(I(J1),J1=1,9),MTYPE  ! Read in element information
        
        !       Save element information
        DO J1=1,9
           XY(2*J1-1,N)=X(I(J1))  ! Coordinates of the element's j node
           XY(2*J1,N)=Y(I(J1))
        END DO
        
        MATP(N)=MTYPE  ! Material type
        
        DO L=1,18
           LM(L,N)=0
        END DO
        
        DO L=1,2
           DO J1=1,9
              LM(L+2*J1-2,N)=ID(L,I(J1))     ! Connectivity matrix
           END DO
        END DO
        
=======
     N=0
     DO WHILE (N .NE. NUME)
        READ (IIN,'(12I5)') N,(I(J1),J1=1,9),MTYPE  ! Read in element information
        
        !       Save element information
        DO J1=1,9
           XY(2*J1-1,N)=X(I(J1))  ! Coordinates of the element's j node
           XY(2*J1,N)=Y(I(J1))
        END DO
        
        MATP(N)=MTYPE  ! Material type
        
        DO L=1,18
           LM(L,N)=0
        END DO
        
        DO L=1,2
           DO J1=1,9
              LM(L+2*J1-2,N)=ID(L,I(J1))     ! Connectivity matrix
           END DO
        END DO
        
>>>>>>> master
        !       Update column heights and bandwidth
        CALL COLHT (MHT,ND,LM(1,N))   

        WRITE (IOUT,"(I5,2X,9(4X,I5),7X,I5)") N,(I(J1),J1=1,9),MTYPE

     END DO

<<<<<<< HEAD
     RETURN
        
! Assemble stucture stiffness matrix
  ELSE IF (IND .EQ. 2) THEN

     !CALL GAUSS(NGP, GAUSSP, WEIGHT)
     GAUSSP=GAUSS(1:NGP,NGP)
     WEIGHT=GAUSS(1:NGP,NGP)
      
     DO N=1,NUME
        MTYPE=MATP(N)
        
        POISSONE=POISSON(MTYPE)
        IF (NPAR(4)==0) THEN
           D(1,1)=1.D0
           D(1,2)=POISSONE
           D(2,1)=POISSONE
           D(2,2)=1.D0
           D(3,3)=(1.D0-POISSONE)/2.D0
           D=E(MTYPE)/(1.D0-POISSONE*POISSONE)*D
        ELSEIF (NPAR(4)==1) THEN
           D(1,1)=1.D0-POISSONE
           D(1,2)=POISSONE
           D(2,1)=POISSONE
           D(2,2)=1.D0-POISSONE
           D(3,3)=(1.D0-2.D0*POISSONE)/2.D0
           D=E(MTYPE)/(1.D0+POISSONE)/(1.D0-2.D0*POISSONE)*D
        END IF
        
        S=0
        DO J1=1,NGP
           XI=GAUSSP(J1)
           DO J2=1,NGP
              ETA=GAUSSP(J2)
              TWOWEIGHT=WEIGHT(J1)*WEIGHT(J2)
              
              NXI(1)=-XI*(1.D0-XI)/2.D0
              NXI(2)=XI*(1.D0+XI)/2.D0
              NXI(3)=1.D0-XI*XI
              NETA(1)=-ETA*(1.D0-ETA)/2.D0
              NETA(2)=ETA*(1.D0+ETA)/2.D0
              NETA(3)=1.D0-ETA*ETA
              
              DNXI(1)=-0.5D0+XI
              DNXI(2)=0.5D0+XI
              DNXI(3)=-2.D0*XI
              DNETA(1)=-0.5D0+ETA
              DNETA(2)=0.5D0+ETA
              DNETA(3)=-2.D0*ETA
              
              JE=0
              DO J3=1,3
                 DO J4=1,3
                    JE(1,1)=JE(1,1)+DNXI(J3)*NETA(J4)*XY(2*POSITION(J3,J4)-1,N)
                    JE(1,2)=JE(1,2)+DNXI(J3)*NETA(J4)*XY(2*POSITION(J3,J4),N)
                    JE(2,1)=JE(2,1)+NXI(J3)*DNETA(J4)*XY(2*POSITION(J3,J4)-1,N)
                    JE(2,2)=JE(2,2)+NXI(J3)*DNETA(J4)*XY(2*POSITION(J3,J4),N)
                 END DO
              END DO
              
              DETJE=JE(1,1)*JE(2,2)-JE(1,2)*JE(2,1)
              
              JENI(1,1)=JE(2,2)/DETJE
              JENI(1,2)=-JE(1,2)/DETJE
              JENI(2,1)=-JE(2,1)/DETJE
              JENI(2,2)=JE(1,1)/DETJE
              
              !DERIVAN(2,9)
              DO J3=1,3
                 DO J4=1,3
                    DERIVAN(1,POSITION(J3,J4))=JENI(1,1)*DNXI(J3)*NETA(J4)+JENI(1,2)*NXI(J3)*DNETA(J4)
                    DERIVAN(2,POSITION(J3,J4))=JENI(2,1)*DNXI(J3)*NETA(J4)+JENI(2,2)*NXI(J3)*DNETA(J4)
                 END DO
              END DO
              
              DO J4=1,9
                 DO J3=1,J4
                    S(2*J3-1,2*J4-1)=S(2*J3-1,2*J4-1)+TWOWEIGHT*DETJE*(DERIVAN(1,J3)*DERIVAN(1,J4)*D(1,1)+DERIVAN(1,J3)*DERIVAN(2,J4)*D(1,3)+DERIVAN(2,J3)*DERIVAN(1,J4)*D(3,1)+DERIVAN(2,J3)*DERIVAN(2,J4)*D(3,3))
                    S(2*J3-1,2*J4)=S(2*J3-1,2*J4)+TWOWEIGHT*DETJE*(DERIVAN(1,J3)*DERIVAN(1,J4)*D(1,3)+DERIVAN(1,J3)*DERIVAN(2,J4)*D(1,2)+DERIVAN(2,J3)*DERIVAN(1,J4)*D(3,3)+DERIVAN(2,J3)*DERIVAN(2,J4)*D(3,2))
                    S(2*J3,2*J4-1)=S(2*J3,2*J4-1)+TWOWEIGHT*DETJE*(DERIVAN(1,J3)*DERIVAN(1,J4)*D(3,1)+DERIVAN(1,J3)*DERIVAN(2,J4)*D(3,3)+DERIVAN(2,J3)*DERIVAN(1,J4)*D(2,1)+DERIVAN(2,J3)*DERIVAN(2,J4)*D(2,3))
                    S(2*J3,2*J4)=S(2*J3,2*J4)+TWOWEIGHT*DETJE*(DERIVAN(1,J3)*DERIVAN(1,J4)*D(3,3)+DERIVAN(1,J3)*DERIVAN(2,J4)*D(3,2)+DERIVAN(2,J3)*DERIVAN(1,J4)*D(2,3)+DERIVAN(2,J3)*DERIVAN(2,J4)*D(2,2))
                 END DO
              END DO
              
           END DO
        END DO
        
        CALL ADDBAN (DA(NP(3)),IA(NP(2)),S,LM(1,N),ND)

     END DO

     RETURN
     
! Stress calculations
  ELSE IF (IND .EQ. 3) THEN
      
=======
     RETURN
        
! Assemble stucture stiffness matrix
  ELSE IF (IND .EQ. 2) THEN

     !CALL GAUSS(NGP, GAUSSP, WEIGHT)
     GAUSSP=GAUSS(1:NGP,NGP)
     WEIGHT=GAUSS(1:NGP,NGP)
      
     DO N=1,NUME
        MTYPE=MATP(N)
        
        POISSONE=POISSON(MTYPE)
        IF (NPAR(4)==0) THEN
           D(1,1)=1.D0
           D(1,2)=POISSONE
           D(2,1)=POISSONE
           D(2,2)=1.D0
           D(3,3)=(1.D0-POISSONE)/2.D0
           D=E(MTYPE)/(1.D0-POISSONE*POISSONE)*D
        ELSEIF (NPAR(4)==1) THEN
           D(1,1)=1.D0-POISSONE
           D(1,2)=POISSONE
           D(2,1)=POISSONE
           D(2,2)=1.D0-POISSONE
           D(3,3)=(1.D0-2.D0*POISSONE)/2.D0
           D=E(MTYPE)/(1.D0+POISSONE)/(1.D0-2.D0*POISSONE)*D
        END IF
        
        S=0
        DO J1=1,NGP
           XI=GAUSSP(J1)
           DO J2=1,NGP
              ETA=GAUSSP(J2)
              TWOWEIGHT=WEIGHT(J1)*WEIGHT(J2)
              
              NXI(1)=-XI*(1.D0-XI)/2.D0
              NXI(2)=XI*(1.D0+XI)/2.D0
              NXI(3)=1.D0-XI*XI
              NETA(1)=-ETA*(1.D0-ETA)/2.D0
              NETA(2)=ETA*(1.D0+ETA)/2.D0
              NETA(3)=1.D0-ETA*ETA
              
              DNXI(1)=-0.5D0+XI
              DNXI(2)=0.5D0+XI
              DNXI(3)=-2.D0*XI
              DNETA(1)=-0.5D0+ETA
              DNETA(2)=0.5D0+ETA
              DNETA(3)=-2.D0*ETA
              
              JE=0
              DO J3=1,3
                 DO J4=1,3
                    JE(1,1)=JE(1,1)+DNXI(J3)*NETA(J4)*XY(2*POSITION(J3,J4)-1,N)
                    JE(1,2)=JE(1,2)+DNXI(J3)*NETA(J4)*XY(2*POSITION(J3,J4),N)
                    JE(2,1)=JE(2,1)+NXI(J3)*DNETA(J4)*XY(2*POSITION(J3,J4)-1,N)
                    JE(2,2)=JE(2,2)+NXI(J3)*DNETA(J4)*XY(2*POSITION(J3,J4),N)
                 END DO
              END DO
              
              DETJE=JE(1,1)*JE(2,2)-JE(1,2)*JE(2,1)
              
              JENI(1,1)=JE(2,2)/DETJE
              JENI(1,2)=-JE(1,2)/DETJE
              JENI(2,1)=-JE(2,1)/DETJE
              JENI(2,2)=JE(1,1)/DETJE
              
              !DERIVAN(2,9)
              DO J3=1,3
                 DO J4=1,3
                    DERIVAN(1,POSITION(J3,J4))=JENI(1,1)*DNXI(J3)*NETA(J4)+JENI(1,2)*NXI(J3)*DNETA(J4)
                    DERIVAN(2,POSITION(J3,J4))=JENI(2,1)*DNXI(J3)*NETA(J4)+JENI(2,2)*NXI(J3)*DNETA(J4)
                 END DO
              END DO
              
              DO J4=1,9
                 DO J3=1,J4
                    S(2*J3-1,2*J4-1)=S(2*J3-1,2*J4-1)+TWOWEIGHT*DETJE*(DERIVAN(1,J3)*DERIVAN(1,J4)*D(1,1)+DERIVAN(1,J3)*DERIVAN(2,J4)*D(1,3)+DERIVAN(2,J3)*DERIVAN(1,J4)*D(3,1)+DERIVAN(2,J3)*DERIVAN(2,J4)*D(3,3))
                    S(2*J3-1,2*J4)=S(2*J3-1,2*J4)+TWOWEIGHT*DETJE*(DERIVAN(1,J3)*DERIVAN(1,J4)*D(1,3)+DERIVAN(1,J3)*DERIVAN(2,J4)*D(1,2)+DERIVAN(2,J3)*DERIVAN(1,J4)*D(3,3)+DERIVAN(2,J3)*DERIVAN(2,J4)*D(3,2))
                    S(2*J3,2*J4-1)=S(2*J3,2*J4-1)+TWOWEIGHT*DETJE*(DERIVAN(1,J3)*DERIVAN(1,J4)*D(3,1)+DERIVAN(1,J3)*DERIVAN(2,J4)*D(3,3)+DERIVAN(2,J3)*DERIVAN(1,J4)*D(2,1)+DERIVAN(2,J3)*DERIVAN(2,J4)*D(2,3))
                    S(2*J3,2*J4)=S(2*J3,2*J4)+TWOWEIGHT*DETJE*(DERIVAN(1,J3)*DERIVAN(1,J4)*D(3,3)+DERIVAN(1,J3)*DERIVAN(2,J4)*D(3,2)+DERIVAN(2,J3)*DERIVAN(1,J4)*D(2,3)+DERIVAN(2,J3)*DERIVAN(2,J4)*D(2,2))
                 END DO
              END DO
              
           END DO
        END DO
        
        CALL ADDBAN (DA(NP(3)),IA(NP(2)),S,LM(1,N),ND)

     END DO

     RETURN
     
! Stress calculations
  ELSE IF (IND .EQ. 3) THEN
      
>>>>>>> master
     GAUSSP=GAUSS(1:NGP,NGP)

     IPRINT=0
     DO N=1,NUME
        IPRINT=IPRINT + 1
        IF (IPRINT.GT.50) IPRINT=1
        IF (IPRINT.EQ.1) WRITE (IOUT,"(//,' S T R E S S  C A L C U L A T I O N S  F O R  ',  &
                                           'E L E M E N T  G R O U P',I4,//,   &
<<<<<<< HEAD
                                           '  ELEMENT',12X,'STRAIN',12X,'STRESS',/,'  NUMBER')") NG 
        MTYPE=MATP(N)
        
        POISSONE=POISSON(MTYPE)
        IF (NPAR(4)==0) THEN
           D(1,1)=1.D0
           D(1,2)=POISSONE
           D(2,1)=POISSONE
           D(2,2)=1.D0
           D(3,3)=(1.D0-POISSONE)/2.D0
           D=E(MTYPE)/(1.D0-POISSONE*POISSONE)*D
        ELSEIF (NPAR(4)==1) THEN
           D(1,1)=1.D0-POISSONE
           D(1,2)=POISSONE
           D(2,1)=POISSONE
           D(2,2)=1.D0-POISSONE
           D(3,3)=(1.D0-2.D0*POISSONE)/2.D0
           D=E(MTYPE)/(1.D0+POISSONE)/(1.D0-2.D0*POISSONE)*D
        END IF
        
        DO J1=1,NGP
           XI=GAUSSP(J1)
           DO J2=1,NGP
              ETA=GAUSSP(J2)
              
              NXI(1)=-XI*(1.D0-XI)/2.D0
              NXI(2)=XI*(1.D0+XI)/2.D0
              NXI(3)=1.D0-XI*XI
              NETA(1)=-ETA*(1.D0-ETA)/2.D0
              NETA(2)=ETA*(1.D0+ETA)/2.D0
              NETA(3)=1.D0-ETA*ETA
              
              DNXI(1)=-0.5D0+XI
              DNXI(2)=0.5D0+XI
              DNXI(3)=-2.D0*XI
              DNETA(1)=-0.5D0+ETA
              DNETA(2)=0.5D0+ETA
              DNETA(3)=-2.D0*ETA
              
              JE=0
              DO J3=1,3
                 DO J4=1,3
                    JE(1,1)=JE(1,1)+DNXI(J3)*NETA(J4)*XY(2*POSITION(J3,J4)-1,N)
                    JE(1,2)=JE(1,2)+DNXI(J3)*NETA(J4)*XY(2*POSITION(J3,J4),N)
                    JE(2,1)=JE(1,1)+NXI(J3)*DNETA(J4)*XY(2*POSITION(J3,J4)-1,N)
                    JE(2,2)=JE(1,1)+NXI(J3)*DNETA(J4)*XY(2*POSITION(J3,J4),N)
                 END DO
              END DO
              
              DETJE=JE(1,1)*JE(2,2)-JE(1,2)*JE(2,1)
              
              JENI(1,1)=JE(2,2)/DETJE
              JENI(1,2)=-JE(1,2)/DETJE
              JENI(2,1)=-JE(2,1)/DETJE
              JENI(2,2)=JE(1,1)/DETJE
              
              !DERIVAN(2,9)
              DO J3=1,3
                 DO J4=1,3
                    DERIVAN(1,POSITION(J3,J4))=JENI(1,1)*DNXI(J3)*NETA(J4)+JENI(1,2)*NXI(J3)*DNETA(J4)
                    DERIVAN(2,POSITION(J3,J4))=JENI(2,1)*DNXI(J3)*NETA(J4)+JENI(2,2)*NXI(J3)*DNETA(J4)
                 END DO
              END DO
              
              STR=0.0
              DO J3=1,9
                 J4=LM(2*J3-1,N)
                 IF (J4.GT.0) THEN
                    STR(1)=STR(1)+DERIVAN(1,J3)*U(J4)
                    STR(3)=STR(3)+DERIVAN(2,J3)*U(J4)
                 END IF
                 
                 J4=LM(2*J3,N)
                 IF (J4.GT.0) THEN
                    STR(2)=STR(2)+DERIVAN(2,J3)*U(J4)
                    STR(3)=STR(3)+DERIVAN(1,J3)*U(J4)
                 END IF
              END DO
              
              P=0.0
              DO J3=1,3
                 DO J4=1,3
                    P(J3)=P(J3)+D(J3,J4)*STR(J4)
                 END DO
              END DO
              
              WRITE (IOUT,"('GAUSS POINT (',E13.6,','E13.6,')')") XI,ETA
              DO J3=1,3
                 WRITE (IOUT,"(17X,E13.6,5X,E13.6)") P(J3),STR(J3)
              END DO
              
           END DO
        END DO
        
=======
                                           '  ELEMENT',12X,'STRAIN',12X,'STRESS',/,'  NUMBER')") NG 
        MTYPE=MATP(N)
        
        POISSONE=POISSON(MTYPE)
        IF (NPAR(4)==0) THEN
           D(1,1)=1.D0
           D(1,2)=POISSONE
           D(2,1)=POISSONE
           D(2,2)=1.D0
           D(3,3)=(1.D0-POISSONE)/2.D0
           D=E(MTYPE)/(1.D0-POISSONE*POISSONE)*D
        ELSEIF (NPAR(4)==1) THEN
           D(1,1)=1.D0-POISSONE
           D(1,2)=POISSONE
           D(2,1)=POISSONE
           D(2,2)=1.D0-POISSONE
           D(3,3)=(1.D0-2.D0*POISSONE)/2.D0
           D=E(MTYPE)/(1.D0+POISSONE)/(1.D0-2.D0*POISSONE)*D
        END IF
        
        DO J1=1,NGP
           XI=GAUSSP(J1)
           DO J2=1,NGP
              ETA=GAUSSP(J2)
              
              NXI(1)=-XI*(1.D0-XI)/2.D0
              NXI(2)=XI*(1.D0+XI)/2.D0
              NXI(3)=1.D0-XI*XI
              NETA(1)=-ETA*(1.D0-ETA)/2.D0
              NETA(2)=ETA*(1.D0+ETA)/2.D0
              NETA(3)=1.D0-ETA*ETA
              
              DNXI(1)=-0.5D0+XI
              DNXI(2)=0.5D0+XI
              DNXI(3)=-2.D0*XI
              DNETA(1)=-0.5D0+ETA
              DNETA(2)=0.5D0+ETA
              DNETA(3)=-2.D0*ETA
              
              JE=0
              DO J3=1,3
                 DO J4=1,3
                    JE(1,1)=JE(1,1)+DNXI(J3)*NETA(J4)*XY(2*POSITION(J3,J4)-1,N)
                    JE(1,2)=JE(1,2)+DNXI(J3)*NETA(J4)*XY(2*POSITION(J3,J4),N)
                    JE(2,1)=JE(1,1)+NXI(J3)*DNETA(J4)*XY(2*POSITION(J3,J4)-1,N)
                    JE(2,2)=JE(1,1)+NXI(J3)*DNETA(J4)*XY(2*POSITION(J3,J4),N)
                 END DO
              END DO
              
              DETJE=JE(1,1)*JE(2,2)-JE(1,2)*JE(2,1)
              
              JENI(1,1)=JE(2,2)/DETJE
              JENI(1,2)=-JE(1,2)/DETJE
              JENI(2,1)=-JE(2,1)/DETJE
              JENI(2,2)=JE(1,1)/DETJE
              
              !DERIVAN(2,9)
              DO J3=1,3
                 DO J4=1,3
                    DERIVAN(1,POSITION(J3,J4))=JENI(1,1)*DNXI(J3)*NETA(J4)+JENI(1,2)*NXI(J3)*DNETA(J4)
                    DERIVAN(2,POSITION(J3,J4))=JENI(2,1)*DNXI(J3)*NETA(J4)+JENI(2,2)*NXI(J3)*DNETA(J4)
                 END DO
              END DO
              
              STR=0.0
              DO J3=1,9
                 J4=LM(2*J3-1,N)
                 IF (J4.GT.0) THEN
                    STR(1)=STR(1)+DERIVAN(1,J3)*U(J4)
                    STR(3)=STR(3)+DERIVAN(2,J3)*U(J4)
                 END IF
                 
                 J4=LM(2*J3,N)
                 IF (J4.GT.0) THEN
                    STR(2)=STR(2)+DERIVAN(2,J3)*U(J4)
                    STR(3)=STR(3)+DERIVAN(1,J3)*U(J4)
                 END IF
              END DO
              
              P=0.0
              DO J3=1,3
                 DO J4=1,3
                    P(J3)=P(J3)+D(J3,J4)*STR(J4)
                 END DO
              END DO
              
              WRITE (IOUT,"('GAUSS POINT (',E13.6,','E13.6,')')") XI,ETA
              DO J3=1,3
                 WRITE (IOUT,"(17X,E13.6,5X,E13.6)") P(J3),STR(J3)
              END DO
              
           END DO
        END DO
        
>>>>>>> master
     END DO

  ELSE 
     STOP " ERROR"
<<<<<<< HEAD
  END IF
        
END SUBROUTINE INE_QUAD
        
=======
  END IF
        
END SUBROUTINE INE_QUAD
        
>>>>>>> master
        