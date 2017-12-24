SUBROUTINE HEXA8
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   To set up storage and call the quadrangle element subroutine    .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    USE GLOBALS
    USE memAllocate
    
    IMPLICIT NONE
    INTEGER :: NUME, NUMMAT, MM, N101, N102, N103, N104, N105, N106, N107

    NUME = NPAR(2)
    NUMMAT = NPAR(3)

! Allocate storage for element group data
  IF (IND == 1) THEN
      MM = 3*NUMMAT*ITWO + 25*NUME + 24*NUME*ITWO
      CALL MEMALLOC(11,"ELEGP",MM,1)
  END IF

  NFIRST=NP(11)   ! Pointer to the first entry in the element group data array
                  ! in the unit of single precision (corresponding to A)

! Calculate the pointer to the arrays in the element group data
! N101: E(NUMMAT)
! N102: POISSON'S RATIO(NUMMAT)
! N103: DEN(NUMMAT) Density of each material
! N104: LM(24,NUME)
! N105: XYZ(24,NUME)
! N106: MTAP(NUME)
  N101=NFIRST
  N102=N101+NUMMAT*ITWO
  N103=N102+NUMMAT*ITWO
  N104=N103+NUMMAT*ITWO
  N105=N104+24*NUME
  N106=N105+24*NUME*ITWO
  N107=N106+NUME
  NLAST=N107

  MIDEST=NLAST - NFIRST

  CALL EXA (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)),DA(NP(4)),IA(NP(5)),   &
       A(N101),A(N102),A(N103),A(N104),A(N105),A(N106))

  RETURN
    
    END SUBROUTINE HEXA8
    
    SUBROUTINE EXA (ID,X,Y,Z,U,MHT,E,POISSON,DEN,LM,XYZ,MATP)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   TRUSS element subroutine                                        .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: ID(6,NUMNP),LM(24,NPAR(2)),MATP(NPAR(2)),MHT(NEQ)
  REAL(8) :: X(NUMNP),Y(NUMNP),Z(NUMNP),E(NPAR(3)),POISSON(NPAR(3)),DEN(NPAR(3)),  &
             XYZ(24,NPAR(2)),U(NEQ)
  REAL(8) :: S(24,24),BE(6,24),D(6,6),NE(8),DISP(24,1),COORD(3,8),ERR!,TEMP(24,24),TEST(24,1)

  INTEGER :: NPAR1, NUME, NUMMAT, ND, I1, I2, I3, I4, I5, I6, I7, I8, L, N
  INTEGER :: MTYPE, IPRINT, I, J, K, SEQ(8)
  REAL(8) :: XX, DETJ, STR(6,8)
  !REAL(8) :: XL2, XL, SQRT, XX, YY, STR, P, V

  NPAR1  = NPAR(1)
  NUME   = NPAR(2)
  NUMMAT = NPAR(3) 

  ND=24

! Read and generate element information
  IF (IND .EQ. 1) THEN

     WRITE (IOUT,"(' E L E M E N T   D E F I N I T I O N',//,  &
                   ' ELEMENT TYPE ',13(' .'),'( NPAR(1) ) . . =',I5,/,   &
                   '     EQ.1, TRUSS ELEMENTS',/,      &
                   '     EQ.2, 4-NODE QUADRILATERAL ELEMENTS',/,  &
                   '     EQ.3, NOT AVAILABLE',/         &
                   '     EQ.4, 8-NODE HEXAHEDRON ELEMENTS',//,      &
                   ' NUMBER OF ELEMENTS.',10(' .'),'( NPAR(2) ) . . =',I5,/)") NPAR1,NUME
     
     IF (NUMMAT.EQ.0) NUMMAT=1

     WRITE (IOUT,"(' M A T E R I A L   D E F I N I T I O N',//,  &
                   ' NUMBER OF DIFFERENT SETS OF MATERIAL',/,  &
                   ' AND CROSS-SECTIONAL  CONSTANTS ',         &
                   4 (' .'),'( NPAR(3) ) . . =',I5,/)") NUMMAT

     WRITE (IOUT,"('  SET       YOUNG''S     POISSON''S       DENSITY',/,  &
                   ' NUMBER     MODULUS      RATIO',/,  &
                   15 X,'E',11X,'NU',11X,'RHO')")

     DO I=1,NUMMAT
        READ (IIN,'(I5,3F10.0)') N,E(N),POISSON(N),DEN(N)  ! Read material information
        WRITE (IOUT,"(I5,4X,E12.5,2X,E14.6,4X,E14.6)") N,E(N),POISSON(N),DEN(N)
     END DO

     WRITE (IOUT,"(//,' E L E M E N T   I N F O R M A T I O N',//,  &
                      ' ELEMENT     NODE     NODE     NODE     NODE     NODE     NODE     NODE     NODE       MATERIAL',/,   &
                      ' NUMBER-N     I1       I2       I3       I4       I5       I6       I7       I8       SET NUMBER')")

     N=0
     DO WHILE (N .NE. NUME)
        READ (IIN,'(11I5)') N,I1,I2,I3,I4,I5,I6,I7,I8,MTYPE  ! Read in element information
!CAUTION:   Automatic generation is possible here using KG
!       Save element information
        XYZ(1,N)=X(I1)  ! Coordinates of the element's first node
        XYZ(2,N)=Y(I1)
        XYZ(3,N)=Z(I1)
        
        XYZ(4,N)=X(I2)  ! Coordinates of the element's second node
        XYZ(5,N)=Y(I2)
        XYZ(6,N)=Z(I2)
        
        XYZ(7,N)=X(I3)  ! Coordinates of the element's third node
        XYZ(8,N)=Y(I3)
        XYZ(9,N)=Z(I3)
        
        XYZ(10,N)=X(I4)  ! Coordinates of the element's fourth node
        XYZ(11,N)=Y(I4)
        XYZ(12,N)=Z(I4)
        
        XYZ(13,N)=X(I5)  ! Coordinates of the element's first node
        XYZ(14,N)=Y(I5)
        XYZ(15,N)=Z(I5)
        
        XYZ(16,N)=X(I6)  ! Coordinates of the element's second node
        XYZ(17,N)=Y(I6)
        XYZ(18,N)=Z(I6)
        
        XYZ(19,N)=X(I7)  ! Coordinates of the element's third node
        XYZ(20,N)=Y(I7)
        XYZ(21,N)=Z(I7)
        
        XYZ(22,N)=X(I8)  ! Coordinates of the element's fourth node
        XYZ(23,N)=Y(I8)
        XYZ(24,N)=Z(I8)
        
        MATP(N)=MTYPE  ! Material type

        DO L=1,24
           LM(L,N)=0
        END DO

        DO L=1,3
           LM(L,N)=ID(L,I1)     ! Connectivity matrix
           LM(L+3,N)=ID(L,I2)
           LM(L+6,N)=ID(L,I3)
           LM(L+9,N)=ID(L,I4)
           LM(L+12,N)=ID(L,I5)
           LM(L+15,N)=ID(L,I6)
           LM(L+18,N)=ID(L,I7)
           LM(L+21,N)=ID(L,I8)
        END DO

!       Update column heights and bandwidth
        CALL COLHT (MHT,ND,LM(1,N))   

        WRITE (IOUT,"(I5,6X,7(I5,4X),I5,7X,I5)") N,I1,I2,I3,I4,I5,I6,I7,I8,MTYPE

!       Generate the gravity

        
        
        
        
        
        
        
        
        
        
        
        !V=XL*DEN(MATP(N))*AREA(MATP(N))
        !
        !CALL GRAVITY(LM(:,N),V)
        
     END DO

     RETURN

! Assemble stucture stiffness matrix
  ELSE IF (IND .EQ. 2) THEN

     DO N=1,NUME
        MTYPE=MATP(N)
        
        ! Generate the D matrix
        D = 0.
            XX=E(MTYPE)/((1.D0+POISSON(MTYPE))*(1.D0-2.D0*POISSON(MTYPE)))       ! XX gives the coefficient for plane strain matrix
            D(1,1) = 1.D0 - POISSON(MTYPE)
            D(2,2) = 1.D0 - POISSON(MTYPE)
            D(3,3) = 1.D0 - POISSON(MTYPE)
            D(1,2) = POISSON(MTYPE)
            D(2,1) = POISSON(MTYPE)
            D(2,3) = POISSON(MTYPE)
            D(3,2) = POISSON(MTYPE)
            D(3,1) = POISSON(MTYPE)
            D(1,3) = POISSON(MTYPE)
            D(4,4) = (1.D0-2.D0*POISSON(MTYPE))/2.D0
            D(5,5) = (1.D0-2.D0*POISSON(MTYPE))/2.D0
            D(6,6) = (1.D0-2.D0*POISSON(MTYPE))/2.D0
        
        S = 0.
        DO I = 1,2
            DO J = 1,2
                DO K = 1,2
                    CALL BEJE_HEXA8(GAUSS(I,2),GAUSS(J,2),GAUSS(K,2), BE, DETJ, XYZ(:,N))
                    !TEMP = MATMUL(TRANSPOSE(BE),MATMUL(D,BE))
                    S = S + GAUSSW(I,2)*GAUSSW(J,2)*GAUSSW(K,2)*DETJ * MATMUL(TRANSPOSE(BE),MATMUL(D,BE))
                END DO
            END DO
        END DO
        
        S = S*XX

        !TEST = 0.D0
        !TEST(1:3:24,1) = 1.D-3
        !TEST = MATMUL(S,TEST)
        
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
                                           '       ELEMENT',23X,'COORDIANTE STRESS IN EACH NODE IN PARTICULAR DIRECTION',/, &
                                           '  NUMBER & DIRECTION',23X,'1',33X,'2',33X,'3',33X,'4',33X,'5',33X,'6',33X,'7',33X,'8')") NG
        MTYPE=MATP(N)

        ! Generate the D matrix
        D = 0.
            XX=E(MTYPE)/((1.D0+POISSON(MTYPE))*(1.D0-2.D0*POISSON(MTYPE)))       ! XX gives the coefficient for plane strain matrix
            D(1,1) = 1.D0 - POISSON(MTYPE)
            D(2,2) = 1.D0 - POISSON(MTYPE)
            D(3,3) = 1.D0 - POISSON(MTYPE)
            D(1,2) = POISSON(MTYPE)
            D(2,1) = POISSON(MTYPE)
            D(2,3) = POISSON(MTYPE)
            D(3,2) = POISSON(MTYPE)
            D(3,1) = POISSON(MTYPE)
            D(1,3) = POISSON(MTYPE)
            D(4,4) = (1.D0-2.D0*POISSON(MTYPE))/2.D0
            D(5,5) = (1.D0-2.D0*POISSON(MTYPE))/2.D0
            D(6,6) = (1.D0-2.D0*POISSON(MTYPE))/2.D0
            D = D*XX
        
        ! Pull down the displacement of every node
        DISP = 0.D0
        DO I = 1,24
            IF(LM(I,N) .NE. 0) DISP(I,1) = U(LM(I,N))
        END DO
        
        S = 0.
        DO I = 1,2
            DO J = 1,2
                DO K = 1,2
                    CALL BEJE_HEXA8(GAUSS(I,2),GAUSS(J,2),GAUSS(K,2), BE, DETJ, XYZ(:,N))
                    STR(:,4*(K-1)+2*(J-1)+I) = MATMUL(BE,DISP(:,1))
                    CALL NE_HEXA8(GAUSS(I,2),GAUSS(J,2),GAUSS(K,2), NE)
                    COORD(:,(4*(K-1)+2*(J-1)+I)) = &
                            NE(1)*XYZ(1:3,N) + NE(2)*XYZ(4:6,N) + NE(3)*XYZ(7:9,N) + NE(4)*XYZ(10:12,N) + &
                            NE(5)*XYZ(13:15,N) + NE(6)*XYZ(16:18,N) + NE(7)*XYZ(19:21,N) + NE(8)*XYZ(22:24,N)
                END DO
            END DO
        END DO
        
        CALL ERR_HEXA8(ERR,COORD,STR,D)
        
        SEQ = (/1,2,4,3,5,6,8,7/)           ! Sequence of the stress shown
        
        WRITE (IOUT,"(6X,I5,'X',17X,15(E13.6,4X),E13.6)") N,(COORD(1,SEQ(I)),STR(1,SEQ(I)), I=1,8)
        WRITE (IOUT,"(6X,I5,'Y',17X,15(E13.6,4X),E13.6)") N,(COORD(2,SEQ(I)),STR(2,SEQ(I)), I=1,8)
        WRITE (IOUT,"(6X,I5,'Z',17X,15(E13.6,4X),E13.6)") N,(COORD(3,SEQ(I)),STR(3,SEQ(I)), I=1,8)
        WRITE (IOUT,"(6X,I5,'YZSHEAR',28X,7(E13.6,21X),E13.6)") N,STR(4,SEQ)
        WRITE (IOUT,"(6X,I5,'ZXSHEAR',28X,7(E13.6,21X),E13.6)") N,STR(5,SEQ)
        WRITE (IOUT,"(6X,I5,'XYSHEAR',28X,7(E13.6,21X),E13.6)") N,STR(6,SEQ)
        
     END DO

     WRITE (IOUT,"(/,6X,'ENERGY ERROR IS',6X,E13.6)") ERR
     
  ELSE 
     STOP "*** ERROR *** Invalid IND value."
  END IF

    END SUBROUTINE EXA

    
    ! Calculate the shape function of element
    SUBROUTINE NE_HEXA8(PSI, ETA, XI, NE)
    
    IMPLICIT NONE
    
    REAL(8), INTENT(IN)             :: PSI, ETA, XI
    REAL(8), INTENT(OUT)            :: NE(8)           ! Shape function of the quadralic elements
    
    NE(1) = 0.125D0*(1.D0-PSI)*(1.D0-ETA)*(1.D0-XI)
    NE(2) = 0.125D0*(1.D0+PSI)*(1.D0-ETA)*(1.D0-XI)
    NE(3) = 0.125D0*(1.D0+PSI)*(1.D0+ETA)*(1.D0-XI)
    NE(4) = 0.125D0*(1.D0-PSI)*(1.D0+ETA)*(1.D0-XI)
    NE(5) = 0.125D0*(1.D0-PSI)*(1.D0-ETA)*(1.D0+XI)
    NE(6) = 0.125D0*(1.D0+PSI)*(1.D0-ETA)*(1.D0+XI)
    NE(7) = 0.125D0*(1.D0+PSI)*(1.D0+ETA)*(1.D0+XI)
    NE(8) = 0.125D0*(1.D0-PSI)*(1.D0+ETA)*(1.D0+XI)
    
    END SUBROUTINE NE_HEXA8
    
    
    ! Calculate the strain matrix and the Jacobi matrix for the element
    SUBROUTINE BEJE_HEXA8(PSI, ETA, XI, BE, DETJ, XYZ)
    
    IMPLICIT NONE
    
    REAL(8), INTENT(IN)             :: PSI, ETA, XI, XYZ(24)
    REAL(8), INTENT(OUT)            :: BE(6,24), DETJ
    REAL(8)                         :: PN(3,8), JE(3,3), JEI(3,3)        ! Partial of Ne matrix
    
    ! Partial ETA,PSI,XI to Ne matrix
    PN(1,:) = 0.125D0*(/ (ETA-1.D0)*(1.D0-XI), (1.D0-ETA)*(1.D0-XI), (1.D0+ETA)*(1.D0-XI), (-1.D0-ETA)*(1.D0-XI),&
                         (ETA-1.D0)*(1.D0+XI), (1.D0-ETA)*(1.D0+XI), (1.D0+ETA)*(1.D0+XI), (-1.D0-ETA)*(1.D0+XI) /)
    PN(2,:) = 0.125D0*(/ (PSI-1.D0)*(1.D0-XI), (-1.D0-PSI)*(1.D0-XI), (1.D0+PSI)*(1.D0-XI), (1.D0-PSI)*(1.D0-XI),&
                         (PSI-1.D0)*(1.D0+XI), (-1.D0-PSI)*(1.D0+XI), (1.D0+PSI)*(1.D0+XI), (1.D0-PSI)*(1.D0+XI) /)
    PN(3,:) = 0.125D0*(/ -(1.D0-PSI)*(1.D0-ETA), -(1.D0+PSI)*(1.D0-ETA), -(1.D0+PSI)*(1.D0+ETA), -(1.D0-PSI)*(1.D0+ETA),&
                         (1.D0-PSI)*(1.D0-ETA), (1.D0+PSI)*(1.D0-ETA), (1.D0+PSI)*(1.D0+ETA), (1.D0-PSI)*(1.D0+ETA) /)
    
    ! Jacobi matrix and its determinant, J = PN*[X Y]
    JE(1,1) = DOT_PRODUCT(PN(1,:),XYZ(1:24:3))
    JE(2,1) = DOT_PRODUCT(PN(2,:),XYZ(1:24:3))
    JE(3,1) = DOT_PRODUCT(PN(3,:),XYZ(1:24:3))
    JE(1,2) = DOT_PRODUCT(PN(1,:),XYZ(2:24:3))
    JE(2,2) = DOT_PRODUCT(PN(2,:),XYZ(2:24:3))
    JE(3,2) = DOT_PRODUCT(PN(3,:),XYZ(2:24:3))
    JE(1,3) = DOT_PRODUCT(PN(1,:),XYZ(3:24:3))
    JE(2,3) = DOT_PRODUCT(PN(2,:),XYZ(3:24:3))
    JE(3,3) = DOT_PRODUCT(PN(3,:),XYZ(3:24:3))
    
    DETJ = JE(3,3)*JE(2,2)*JE(1,1) + JE(1,2)*JE(2,3)*JE(3,1) + JE(1,3)*JE(2,1)*JE(3,2)&
         - JE(1,3)*JE(2,2)*JE(3,1) - JE(1,2)*JE(2,1)*JE(3,3) - JE(1,1)*JE(2,3)*JE(3,2)
    
    IF(ABS(DETJ) .LE. EPSILON(DETJ)) THEN
        PRINT *, '***ERROR***  ELEMENT JACOBI MATRIX SINGULAR'
        READ (*,*)
    END IF
        
    JEI(1,1) = JE(2,2)*JE(3,3) - JE(2,3)*JE(3,2)
    JEI(2,2) = JE(1,1)*JE(3,3) - JE(1,3)*JE(3,1)
    JEI(3,3) = JE(1,1)*JE(2,2) - JE(1,2)*JE(2,1)
    JEI(1,2) = -(JE(1,2)*JE(3,3)-JE(3,2)*JE(1,3))
    JEI(2,1) = -(JE(2,1)*JE(3,3)-JE(3,1)*JE(2,3))
    JEI(1,3) = (JE(1,2)*JE(2,3)-JE(2,2)*JE(1,3))
    JEI(3,1) = (JE(2,1)*JE(3,2)-JE(2,2)*JE(3,1))
    JEI(2,3) = -(JE(1,1)*JE(2,3)-JE(2,1)*JE(1,3))
    JEI(3,2) = -(JE(1,1)*JE(3,2)-JE(3,1)*JE(1,2))
    JEI = JEI/DETJ               ! Inverse matrix
    
    ! Partial x,y,z to Ne matrix, = (J^-1)*PN
    PN = MATMUL(JEI,PN)
    
    ! BE = PARTIAL X,Y,Z TO NE IN STRAIN WAY, BE(EX,EY,EZ, EYZ, EZX, EXY)
    BE = 0.
    BE(1,1:24:3) = PN(1,:)
    BE(2,2:24:3) = PN(2,:)
    BE(3,3:24:3) = PN(3,:)
    BE(4,3:24:3) = PN(2,:)
    BE(4,2:24:3) = PN(3,:)
    BE(5,1:24:3) = PN(3,:)
    BE(5,3:24:3) = PN(1,:)
    BE(6,1:24:3) = PN(2,:)
    BE(6,2:24:3) = PN(1,:)
    
    END SUBROUTINE BEJE_HEXA8
    
    
    ! Calculate the error for each point
    SUBROUTINE ERR_HEXA8(ERR,COORD,STR,D)
    
    
    USE GLOBALS, ONLY : NPAR
    IMPLICIT NONE
    
    REAL(8), INTENT(IN)     :: STR(6,8), COORD(3,8),D(6,6)
    REAL(8), INTENT(INOUT)  :: ERR
    REAL(8)                 :: EXA(6,8),ERRT
    INTEGER                 :: I
    
    EXA = 0.D0
    EXA(1,:) = 0.002D0*COORD(1,:)**2
    
    EXA = EXA - STR
    DO I = 1,4
        ERRT = DOT_PRODUCT( (EXA(:,I)) , MATMUL(D,EXA(:,I)) )
        ERR = ERR + SQRT(ERRT)/NPAR(2)
    END DO
    
    END SUBROUTINE
