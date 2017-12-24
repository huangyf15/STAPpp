  SUBROUTINE INFINITE
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                 .
! .   To set up storage and call the infinite element subroutine    .
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
      MM = 2*NUMMAT*ITWO + 9*NUME + 8*NUME*ITWO
      CALL MEMALLOC(11,"ELEGP",MM,1)
  END IF

  NFIRST=NP(11)   ! Pointer to the first entry in the element group data array
                  ! in the unit of single precision (corresponding to A)

! Calculate the pointer to the arrays in the element group data
! N101: E(NUMMAT)
! N102: POSSION(NUMMAT)
! N103: LM(8,NUME)
! N104: XY(8,NUME)
! N105: MTAP(NUME)
  N101=NFIRST
  N102=N101+NUMMAT*ITWO
  N103=N102+NUMMAT*ITWO
  N104=N103+8*NUME
  N105=N104+8*NUME*ITWO
  N106=N105+NUME
  NLAST=N106

  MIDEST=NLAST - NFIRST

  CALL SUBINFI (IA(NP(1)),DA(NP(2)),DA(NP(3)),DA(NP(4)),IA(NP(5)),   &
       A(N101),A(N102),A(N103),A(N104),A(N105))

  RETURN

  END SUBROUTINE INFINITE


SUBROUTINE SUBINFI (ID,X,Y,U,MHT,E,POSSION,LM,XY,MATP)
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
! .                                                                   .
! .   infinite element subroutine                                     .
! .                                                                   .
! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  USE GLOBALS
  USE MEMALLOCATE

  IMPLICIT NONE
  INTEGER :: ID(6,NUMNP),LM(8,NPAR(2)),MATP(NPAR(2)),MHT(NEQ)
  REAL(8) :: X(NUMNP),Y(NUMNP),E(NPAR(3)),POSSION(NPAR(3)),  &
             XY(8,NPAR(2)),U(NEQ)
  REAL(8) :: S(8,8),ST(3,3),GN(2,4),GM(2,4),JAC(2,2),INVJ(2,2),BE(3,8), &
             D(8),STRAIN(3),STRESS(3)

  INTEGER :: NPAR1, NUME, NUMMAT, ND, I1, I2, I3, I4
  INTEGER :: I, J, L, N, II, JJ, KK, LL
  INTEGER :: MTYPE, IPRINT
  REAL(8) :: GP, ETA, XI, ABSJ

  NPAR1  = NPAR(1)
  NUME   = NPAR(2)
  NUMMAT = NPAR(3) 

  ND=8

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
                      ' ELEMENT     NODE     NODE     NODE     NODE    MATERIAL',/,   &
                      ' NUMBER-N     I1       I2       I3       I4    SET NUMBER')")

     N=0
     DO WHILE (N .NE. NUME)
        READ (IIN,'(6I5)') N,I1,I2,I3,I4,MTYPE  ! Read in element information

!       Save element information
        XY(1,N)=X(I1)  ! Coordinates of the element's  node
        XY(2,N)=Y(I1)
        XY(3,N)=X(I2)
        XY(4,N)=Y(I2)
        XY(5,N)=X(I3)
        XY(6,N)=Y(I3)
        XY(7,N)=X(I4)
        XY(8,N)=Y(I4)

        MATP(N)=MTYPE  ! Material type

        DO L=1,8
           LM(L,N)=0
        END DO

        DO L=1,2
           LM(L,N)=ID(L,I1)     ! Connectivity matrix
           LM(L+2,N)=ID(L,I2)
           LM(L+4,N)=ID(L,I3)
           LM(L+6,N)=ID(L,I4)
        END DO

!       Update column heights and bandwidth
        CALL COLHT (MHT,ND,LM(1,N))   

        WRITE (IOUT,"(I5,6X,I5,4X,I5,4X,I5,4X,I5,4X,I5)") N,I1,I2,I3,I4,MTYPE

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
        
        GP = 1/SQRT(3.D0)   ! Gauss point
        
        ! Get the element stiffness matrix K_e
        S = 0.
        DO I=1,2
            DO J=1,2
                JAC = 0.
                BE  = 0.
                XI  = (-3+2*I)*GP
                ETA = (-3+2*J)*GP
                
                GN(1,1) = -(ETA+1)/4    ! Ger GN matrix
                GN(1,2) =  (ETA-1)/4
                GN(1,3) = -(ETA-1)/4
                GN(1,4) =  (ETA+1)/4
                GN(2,1) = -(XI-1)/4
                GN(2,2) =  (XI-1)/4
                GN(2,3) = -(XI+1)/4
                GN(2,4) =  (XI+1)/4
                
                GM(1,1) = -(1+ETA)/((1-XI)*(1-XI))    ! Ger GM matrix
                GM(1,2) = -(1-ETA)/((1-XI)*(1-XI))
                GM(1,3) =  (1-ETA)/((1-XI)*(1-XI))
                GM(1,4) =  (1+ETA)/((1-XI)*(1-XI))
                GM(2,1) = -XI/(1-XI)
                GM(2,2) =  XI/(1-XI)
                GM(2,3) = -(1+XI)/(2*(1-XI))
                GM(2,4) =  (1+XI)/(2*(1-XI))
                ! Get Jacobi matrix
                DO II=1,4
                    JAC(1,1) = JAC(1,1)+GM(1,II)*XY(2*II-1,N)
                    JAC(1,2) = JAC(1,2)+GM(1,II)*XY(2*II,N)
                    JAC(2,1) = JAC(2,1)+GM(2,II)*XY(2*II-1,N)
                    JAC(2,2) = JAC(2,2)+GM(2,II)*XY(2*II,N)
                END DO
                ABSJ = JAC(1,1)*JAC(2,2)-JAC(1,2)*JAC(2,1)
                INVJ(1,1) = JAC(2,2)    ! The inverse of Jacobi matrix
                INVJ(1,2) =-JAC(1,2)
                INVJ(2,1) =-JAC(2,1)
                INVJ(2,2) = JAC(1,1)
                ! Get the B_e matrix
                DO II=1,4
                    DO JJ=1,2
                        BE(JJ,2*II+JJ-2) = INVJ(JJ,1)*GN(1,II)+INVJ(JJ,2)*GN(2,II)
                        BE(3,2*II-JJ+1)  = BE(JJ,2*II+JJ-2)
                    END DO
                END DO
                ! Update the K_e matrix
                DO II=1,8
                    DO JJ=1,8
                        DO KK=1,3
                            DO LL=1,3
                                S(II,JJ)=S(II,JJ)+ST(LL,KK)*BE(LL,II)*BE(KK,JJ)/ABSJ
                            END DO
                        END DO
                    END DO
                END DO
            END DO
        END DO
        
        CALL ADDBAN (DA(NP(3)),IA(NP(2)),S,LM(1,N),ND)

     END DO

     RETURN

! Stress calculations
  ELSE IF (IND .EQ. 3) THEN

     IPRINT=0
     DO N=1,NUME
        IPRINT=IPRINT + 1
        IF (IPRINT.GT.50) IPRINT=1
        IF (IPRINT.EQ.1) WRITE (IOUT,"(//,' S T R E S S  C A L C U L A T I O N S  F O R  '//,  &
                                           'E L E M E N T  G R O U P',I4,'  A T  G A U S S P O I N T'//,   &
                                           '  ELEMENT',11X,'SIGMA_XX',9X,'SIGMA_YY',9X,'GAMMA_XY',/,'  NUMBER')") NG
        MTYPE=MATP(N)
        
        DO L=1,8
            I   =LM(L,N)
            IF (I .NE. 0) THEN
                D(L)=U(I)
            ELSE
                D(L)=0
            END IF
        END DO
        
        ! Set the constitutive matrix
        ST(1,1)=E(MTYPE)/(1-POSSION(MTYPE)*POSSION(MTYPE))
        ST(2,2)=ST(1,1)
        ST(1,2)=POSSION(MTYPE)*ST(1,1)
        ST(2,1)=ST(1,2)
        ST(3,3)=ST(1,1)*(1-POSSION(MTYPE))/2
        
        GP = 1/SQRT(3.D0)   ! Gauss point
        
        ! Get the element stiffness matrix K_e
        S = 0.
        DO I=1,2
            DO J=1,2
                JAC = 0.
                BE  = 0.
                XI  = (-3+2*I)*GP
                ETA = (-3+2*J)*GP
                
                GN(1,1) = -(ETA+1)/4    ! Ger GN matrix
                GN(1,2) =  (ETA-1)/4
                GN(1,3) = -(ETA-1)/4
                GN(1,4) =  (ETA+1)/4
                GN(2,1) = -(XI-1)/4
                GN(2,2) =  (XI-1)/4
                GN(2,3) = -(XI+1)/4
                GN(2,4) =  (XI+1)/4
                
                GM(1,1) = -(1+ETA)/((1-XI)*(1-XI))    ! Ger GM matrix
                GM(1,2) = -(1-ETA)/((1-XI)*(1-XI))
                GM(1,3) =  (1-ETA)/((1-XI)*(1-XI))
                GM(1,4) =  (1+ETA)/((1-XI)*(1-XI))
                GM(2,1) = -XI/(1-XI)
                GM(2,2) =  XI/(1-XI)
                GM(2,3) = -(1+XI)/(2*(1-XI))
                GM(2,4) =  (1+XI)/(2*(1-XI))
                ! Get Jacobi matrix
                DO II=1,4
                    JAC(1,1) = JAC(1,1)+GM(1,II)*XY(2*II-1,N)
                    JAC(1,2) = JAC(1,2)+GM(1,II)*XY(2*II,N)
                    JAC(2,1) = JAC(2,1)+GM(2,II)*XY(2*II-1,N)
                    JAC(2,2) = JAC(2,2)+GM(2,II)*XY(2*II,N)
                END DO
                ABSJ = JAC(1,1)*JAC(2,2)-JAC(1,2)*JAC(2,1)
                INVJ(1,1) = JAC(2,2)    ! The inverse of Jacobi matrix
                INVJ(1,2) =-JAC(1,2)
                INVJ(2,1) =-JAC(2,1)
                INVJ(2,2) = JAC(1,1)
                ! Get the B_e matrix
                DO II=1,4
                    DO JJ=1,2
                        BE(JJ,2*II+JJ-2) = INVJ(JJ,1)*GN(1,II)+INVJ(JJ,2)*GN(2,II)
                        BE(3,2*II-JJ+1)  = BE(JJ,2*II+JJ-2)
                    END DO
                END DO
                
                
                ! Calculate the strain in gauss point
                STRAIN = 0.
                STRESS = 0.
                DO KK=1,3
                   DO II=1,8
                       STRAIN(KK)=STRAIN(KK)+BE(KK,II)*D(II)/ABSJ
                   END DO
                END DO
                ! Calculate the stress in gauss point
                DO KK=1,3
                   DO II=1,3
                       STRESS(KK)=STRESS(KK)+ST(KK,II)*STRAIN(II)
                   END DO
                END DO
                WRITE (IOUT,"(1X,I5,11X,E13.6,4X,E13.6,4X,E13.6)") N,STRESS(1),STRESS(2),STRESS(3)
            END DO
        END DO
     END DO

  ELSE 
     STOP "ERROR"
  END IF

END SUBROUTINE SUBINFI
