SUBROUTINE TRIANGLE
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                 .
! .   To set up storage and call the 3T element subroutine          .
! .                                                                 .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: NUME, NUMMAT, MM, N101, N102, N103, N104, N105, N106

  NUME = NPAR(2)
  NUMMAT = NPAR(3)

! Allocate storage for element group data
  IF (IND == 1) THEN
      MM = 2*NUMMAT*ITWO + 7*NUME + 6*NUME*ITWO
      CALL MEMALLOC(11,"ELEGP",MM,1)
  END IF
  N101=NFIRST
  N102=N101+NUMMAT*ITWO
  N103=N102+NUMMAT*ITWO
  N104=N103+6*NUME
  N105=N104+6*NUME*ITWO
  N106=N105+NUME
  NLAST=N106

  MIDEST=NLAST - NFIRST

  CALL SUBTRI (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)),IA(NP(5)),   &
       A(N101),A(N102),A(N103),A(N104),A(N105))

  RETURN

END SUBROUTINE TRIANGLE


SUBROUTINE SUBTRI (ID,X,Y,U,MHT,E,POSSION,LM,XY,MATP)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   3T element subroutine                                           .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: ID(6,NUMNP),LM(6,NPAR(2)),MATP(NPAR(2)),MHT(NEQ)
  REAL(8) :: X(NUMNP),Y(NUMNP),E(NPAR(3)),POSSION(NPAR(3)),  &
             XY(6,NPAR(2)),U(NEQ)
  REAL(8) :: S(6,6),ST(3,3),BE(3,6), D(6),STRAIN(3),STRESS(3)

  INTEGER :: NPAR1, NUME, NUMMAT, ND, I1, I2, I3, I4
  INTEGER :: I, J, L, N, II, JJ, KK, LL
  INTEGER :: MTYPE, IPRINT
  REAL(8) :: AREA

  NPAR1  = NPAR(1)
  NUME   = NPAR(2)
  NUMMAT = NPAR(3) 

  ND=6

! Read and generate element information
  IF (IND .EQ. 1) THEN

     WRITE (IOUT,"(' E L E M E N T   D E F I N I T I O N',//,  &
                   ' ELEMENT TYPE ',13(' .'),'( NPAR(1) ) . . =',I5,/,   &
                   '     EQ.1, 4Q ELEMENTS   ',/,      &
                   '     EQ.2, ELEMENTS CURRENTLY',/,  &
                   '     EQ.3, NOT AVAILABLE',//,      &
                   ' NUMBER OF ELEMENTS.',10(' .'),'( NPAR(2) ) . . =',I5,/)") NPAR1,NUME

     IF (NUMMAT.EQ.0) NUMMAT=1

     WRITE (IOUT,"(' M A T E R I A L   D E F I N I T I O N',//,  &
                   ' NUMBER OF DIFFERENT SETS OF MATERIAL',/,  &
                   ' AND CROSS-SECTIONAL  CONSTANTS ',         &
                   4 (' .'),'( NPAR(3) ) . . =',I5,/)") NUMMAT

     WRITE (IOUT,"('  SET       YOUNG''S         POSSION',/,  &
                   ' NUMBER     MODULUS',10X,'RATIO',/,  &
                   15 X,'E',14X,'V')")

     DO I=1,NUMMAT
        READ (IIN,'(I5,2F10.0)') N,E(N),POSSION(N)  ! Read material information
        WRITE (IOUT,"(I5,4X,E12.5,2X,E14.6)") N,E(N),POSSION(N)
     END DO

     WRITE (IOUT,"(//,' E L E M E N T   I N F O R M A T I O N',//,  &
                      ' ELEMENT     NODE     NODE     NODE    MATERIAL',/,   &
                      ' NUMBER-N     I1       I2       I3    SET NUMBER')")

     N=0
     DO WHILE (N .NE. NUME)
        READ (IIN,'(5I5)') N,I1,I2,I3,MTYPE  ! Read in element information

!       Save element information
        XY(1,N)=X(I1)  ! Coordinates of the element's  node
        XY(2,N)=Y(I1)
        XY(3,N)=X(I2)
        XY(4,N)=Y(I2)
        XY(5,N)=X(I3)
        XY(6,N)=Y(I3)
        MATP(N)=MTYPE  ! Material type

        DO L=1,6
           LM(L,N)=0
        END DO

        DO L=1,2
           LM(L,N)=ID(L,I1)     ! Connectivity matrix
           LM(L+2,N)=ID(L,I2)
           LM(L+4,N)=ID(L,I3)
        END DO

!       Update column heights and bandwidth
        CALL COLHT (MHT,ND,LM(1,N))   

        WRITE (IOUT,"(I5,6X,I5,4X,I5,4X,I5,4X,I5)") N,I1,I2,I3,MTYPE

     END DO

     RETURN

! Assemble stucture stiffness matrix
  ELSE IF (IND .EQ. 2) THEN

     DO N=1,NUME
        MTYPE=MATP(N)

        ! Set the constitutive matrix
        ST = 0.
        ST(1,1)=E(MTYPE)/(1-POSSION(MTYPE)*POSSION(MTYPE))
        ST(2,2)=ST(1,1)
        ST(1,2)=POSSION(MTYPE)*ST(1,1)
        ST(2,1)=ST(1,2)
        ST(3,3)=ST(1,1)*(1-POSSION(MTYPE))/2
        
        ! Get the element stiffness matrix K_e
        S = 0.
        ! Get the B_e matrix
        BE  = 0.
        BE(1,1) = XY(4,N)-XY(6,N)
        BE(1,3) = XY(6,N)-XY(2,N)
        BE(1,5) = XY(2,N)-XY(4,N)
        BE(2,2) = XY(5,N)-XY(3,N)
        BE(2,4) = XY(1,N)-XY(5,N)
        BE(2,6) = XY(3,N)-XY(1,N)
        BE(3,1) = BE(2,2)
        BE(3,2) = BE(1,1)
        BE(3,3) = BE(2,4)
        BE(3,4) = BE(1,3)
        BE(3,5) = BE(2,6)
        BE(3,6) = BE(1,5)
        AREA = XY(1,N)*XY(4,N)+XY(3,N)*XY(6,N)+XY(5,N)*XY(2,N)-XY(1,N)*XY(6,N)-XY(3,N)*XY(2,N)-XY(5,N)*XY(4,N)
        AREA = AREA/2
        ! Update the K_e matrix
        DO I=1,6
            DO J=1,6
                DO II=1,3
                    DO JJ=1,3
                        S(I,J)=S(I,J)+ST(JJ,II)*BE(JJ,I)*BE(II,J)/(4*AREA)
                    ENDDO
                ENDDO
            ENDDO
        ENDDO
        
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
                                           'E L E M E N T  G R O U P',I4,'A T  G A U S S P O I N T'//,   &
                                           '  ELEMENT',11X,'SIGMA_XX',9X,'SIGMA_YY',9X,'GAMMA_XY',/,'  NUMBER')") NG
        MTYPE=MATP(N)
        
        DO L=1,6
            I   =LM(L,N)
            IF (I .NE. 0) THEN
                D(L)=U(I)
            ELSE
                D(L)=0
            ENDIF
        ENDDO
        
        ! Set the constitutive matrix
        ST = 0.
        ST(1,1)=E(MTYPE)/(1-POSSION(MTYPE)*POSSION(MTYPE))
        ST(2,2)=ST(1,1)
        ST(1,2)=POSSION(MTYPE)*ST(1,1)
        ST(2,1)=ST(1,2)
        ST(3,3)=ST(1,1)*(1-POSSION(MTYPE))/2
        
        ! Get the element stiffness matrix K_e
        S = 0.
        ! Get the B_e matrix
        BE  = 0.
        BE(1,1) = XY(4,N)-XY(6,N)
        BE(1,3) = XY(6,N)-XY(2,N)
        BE(1,5) = XY(2,N)-XY(4,N)
        BE(2,2) = XY(5,N)-XY(3,N)
        BE(2,4) = XY(1,N)-XY(5,N)
        BE(2,6) = XY(3,N)-XY(1,N)
        BE(3,1) = BE(2,2)
        BE(3,2) = BE(1,1)
        BE(3,3) = BE(2,4)
        BE(3,4) = BE(1,3)
        BE(3,5) = BE(2,6)
        BE(3,6) = BE(1,5)
        AREA = XY(1,N)*XY(4,N)+XY(3,N)*XY(6,N)+XY(5,N)*XY(2,N)-XY(1,N)*XY(6,N)-XY(3,N)*XY(2,N)-XY(5,N)*XY(4,N)

        ! Calculate the strain in gauss point
        STRAIN = 0.
        STRESS = 0.
        DO I=1,3
            DO II=1,6
                STRAIN(I)=STRAIN(I)+BE(I,II)*D(II)/AREA
            ENDDO
        ENDDO
        ! Calculate the stress in gauss point
        DO L=1,3
            DO II=1,3
                STRESS(L)=STRESS(L)+ST(L,II)*STRAIN(II)
            ENDDO
        ENDDO
        WRITE (IOUT,"(1X,I5,11X,E13.6,4X,E13.6,4X,E13.6)") N,STRESS(1),STRESS(2),STRESS(3)

     END DO

  ELSE 
     STOP "ERROR "
  END IF

END SUBROUTINE SUBTRI
