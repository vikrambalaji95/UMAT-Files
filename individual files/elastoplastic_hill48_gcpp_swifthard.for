      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
	  
      REAL*8 STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)

! ----------------------------------------------------------------
!     UMAT FOR HILL 48 YIELD CRITERION
!     SWIFT HARDENING LAW - Y = K*(EPS0 + EQPLAS)^N
!     IMPLICIT GENERALISED CLOSED POINT PROJECTION ALGORITHM
! ----------------------------------------------------------------

      REAL*8 E, NU, K, EPS0, N_EXP, EQPLAS, LAMBDA, MU, HILL,
     1 ELSTIFF(NTENS, NTENS), ELCOMP(NTENS, NTENS), YIELDFN, HARD(3),
     2 DEPLAS(NTENS), DGAMMA, DDGAMMA, DEN, B(NTENS), BB(NTENS, NTENS),
     3 FLOW(NTENS), RESIDUAL(NTENS), P(NTENS, NTENS), HES(NTENS, NTENS),
     4 INVHES(NTENS, NTENS), TERM1, TERM2, DELSTR(NTENS),
     5 TERM3(NTENS, NTENS), F, G, H, L, M, N, R0, R45, R90
      
      PARAMETER(TOLER = 1.D-6)

!     MATERIAL PROPERTIES INPUT

      ! Elastic Properties
      E = PROPS(1) ! Young's modulus
      NU = PROPS(2) ! Poisson's ratio
      
      ! Swift hardening parameters
      K = PROPS(3)
      EPS0 = PROPS(4)
      N_EXP = PROPS(5) ! Hardening exponent
      
      ! Anisotropy (Lankford) coefficients
      R0 = PROPS(6)
      R45 = PROPS(7)
      R90 = PROPS(8)
    
!     RECOVER SOLUTION DEPENDENT STATE VARIABLES FROM PREVIOUS STEP TIME INCREMENT
      EQPLAS = STATEV(1) ! EQUIVALENT PLASTIC STRAIN

!     CALCULATE HILL YIELD FUNCTION COEFFICIENTS

      F = R0/(R90*(1.0 + R0))
      G = 1.0/(1.0 + R0)
      H = R0/(1.0 + R0)
      N = (R0 + R90)*(2.0*R45 - 1.0)/(2*R90*(1.0 + R0))
      L = 1.5
      M = 1.5
      
!     SET UP P MATRIX
      
      P = 0.0 ! Initialise
      
      IF (NDI.EQ.2) THEN ! Plane stress
          
          P(1,1) = G + H
          P(1,2) = - H
          P(2,1) = -H
          P(2,2) = F + H
          P(3,3) = 2.0*N
          
      ELSE ! 3D/ Plane strain
      
          P(1,1) = G + H
          P(1,2) = -H
          P(1,3) = -G
          P(2,1) = -H
          P(2,2) = F + H
          P(2,3) = - F
          P(3,1) = -G
          P(3,2) = -F
          P(3,3) = G + F
          P(4,4) = 2.0*N
          
          IF (NSHR.EQ.3) THEN
              P(5,5) = 2.0*M
              P(6,6) = 2.0*L
          END IF
      
      END IF
      
!     SET UP ELASTIC STIFFNESS MATRIX
      
      LAMBDA = NU*E/((1.0 + NU)*(1.0 - 2.0*NU)) ! 1ST LAMÉ PARAMETER
      MU = E/(2.0*(1.0 + NU)) ! 2ND LAMÉ PARAMETER

      ELSTIFF = 0.0 ! INITIALISE THE ELASTIC STIFFNESS MATRIX TO ZERO

      IF (NDI.EQ.3) THEN ! NUMBER OF DIRECT STRESS COMPONENTS IS 3
!     THE PROBLEM IS EITHER 3D OR PLANE STRAIN
      DO I = 1, NDI
          DO J = 1, NDI
              ELSTIFF(J, I) = LAMBDA
          END DO
          ELSTIFF(I, I) = LAMBDA + 2.0*MU
      END DO
      DO I = NDI+1, NTENS
          ELSTIFF(I, I) = MU
      END DO

      ELSE ! NUMBER OF DIRECT STRESS COMPONENTS IS 2
!     THE PROBLEM IS PLANE STRESS
      DO I = 1, NDI 
          DO J = 1, NDI
              ELSTIFF(J, I) = NU*E/(1.0 - (NU)**2.0)
          END DO
          ELSTIFF(I, I) = E/(1.0 - (NU)**2.0)
      END DO
      DO I = NDI+1, NTENS
          ELSTIFF(I, I) = MU
      END DO

      END IF

!     CALCULATE PREDICTOR/TRIAL STRESS ASSUMING ELASTIC BEHAVIOUR

      STRESS = STRESS + MATMUL(ELSTIFF, DSTRAN)

!     CALCULATE PREDICTED VALUE OF HILL EQUIVALENT STRESS
      
      CALL INNERPROD(STRESS, MATMUL(P, STRESS), HILL, NTENS)
      HILL = SQRT(HILL)

      CALL UHARD(SYIELD,HARD,EQPLAS,EQPLASRT,TIME,DTIME,TEMP,
     1     DTEMP,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,CMNAME,NSTATV,
     2     STATEV,NUMFIELDV,PREDEF,DPRED,NUMPROPS,PROPS)

      YIELDFN = HILL - SYIELD

      IF (YIELDFN.LT.0.0) THEN
       
!     THE STATE OF STRESS IS ELASTIC
!     PREDICTED VALUE OF STRESS IS EQUAL TO THE TRUE STRESS
!     SET ELASTIC STIFFNESS MATRIX AS MATERIAL TANGENT STIFFNESS MATRIX
      DDSDDE = ELSTIFF

      ELSE

!     PLASTIC DEFORMATION HAS COMMENCED
!     SOLVE FOR EQUIVALENT PLASTIC STRAIN INCREMENT USING THE RETURN MAPPING PROCEDURE

      DGAMMA = 0.0 ! INITIALISE THE INCREMENT IN CONSISTENCY PARAMETER TO ZERO

      DEPLAS = 0.0 ! INITIALISE THE INCREMENT IN PLASTIC STRAIN COMPONENTS TO ZERO
      
      DO WHILE (ABS(YIELDFN).GE.TOLER) ! THE ITERATIONS LAST UNTIL THE ABSOLUTE VALUE OF YIELD FUNCTION BECOMES CLOSE TO ZERO 

      ! FLOW VECTOR FOR HILL 48
      
      B = MATMUL(P, STRESS)

      FLOW = B/HILL

      RESIDUAL = DGAMMA*FLOW - DEPLAS ! PLASTIC FLOW RESIDUAL

!     HESSIAN MATRIX FOR VON MISES

      CALL OUTERPROD(B, B, BB, NTENS, NTENS)
      
      CALL INVERSE(ELSTIFF, ELCOMP, NTENS) ! INVERSE OF ELASTIC STIFFNESS MATRIX AKA ELASTIC COMPLIANCE MATRIX

      INVHES = ELCOMP + (DGAMMA/HILL)*(P - (1.0/HILL**2.0)*BB) ! FORMULA FOR INVERSE HESSIAN MATRIX FOR VON MISES

      CALL INVERSE(INVHES, HES, NTENS) ! INVERSE OF INVERSE HESSIAN MATRIX GIVES HESSIAN MATRIX
 
!     OBTAIN INCREMENT IN CONSISTENCY PARAMETER

      CALL INNERPROD(RESIDUAL, MATMUL(HES, FLOW), TERM1, NTENS)
      CALL INNERPROD(FLOW, MATMUL(HES, FLOW), TERM2, NTENS)

      DDGAMMA = (YIELDFN - TERM1)/(HARD(1) + TERM2)

!     UPDATE THE STRESS VALUES

      DELSTR = MATMUL(HES, RESIDUAL + DDGAMMA*FLOW)

      DEPLAS = DEPLAS + MATMUL(ELCOMP, DELSTR)

      DGAMMA = DGAMMA + DDGAMMA

      STRESS = STRESS - DELSTR

      CALL INNERPROD(STRESS, MATMUL(P, STRESS), HILL, NTENS)
      HILL = SQRT(HILL)
      CALL UHARD(SYIELD,HARD,EQPLAS+DGAMMA,EQPLASRT,TIME,DTIME,TEMP,
     1     DTEMP,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,CMNAME,NSTATV,
     2     STATEV,NUMFIELDV,PREDEF,DPRED,NUMPROPS,PROPS)

      YIELDFN = HILL - SYIELD

      END DO

!     EVALUATE THE MATERIAL TANGENT STIFFNESS MATRIX/ JACOBIAN

      CALL OUTERPROD(MATMUL(HES, FLOW), MATMUL(HES, FLOW), 
     1 TERM3, NTENS, NTENS)

      CALL INNERPROD(FLOW, MATMUL(HES,FLOW), DEN, NTENS)

      DDSDDE = HES - (1.0/(DEN + HARD(1)))*TERM3

      END IF

!     UPDATE THE SOLUTION DEPENDENT STATE VARIABLES
      STATEV(1) = EQPLAS + DGAMMA
      STATEV(2) = HILL

      RETURN
      END SUBROUTINE UMAT
! -------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE INVERSE(MAT, INVMAT, NDIM)
      ! SUBROUTINE TO CALCULATE THE INVERSE OF A NDIM X NDIM SQUARE MATRIX USING GAUSS-JORDAN ELIMINATION METHOD
      
      INCLUDE 'ABA_PARAM.INC'
      
      INTEGER NDIM, BIG

      REAL*8 MAT(NDIM, NDIM), INVMAT(NDIM, NDIM), TEMP(NDIM, NDIM), 
     1 PIVOT, FACTOR, DUMMY

      !INITIALISE THE INVERSE MATRIX INVMAT AS IDENTITY MATRIX
      INVMAT = 0.0
      DO I = 1, NDIM
          INVMAT(I, I) = 1.0
      END DO
      
      !STORING OF MAT(NDIM, NDIM) IN TEMP(NDIM, NDIM) TO AVOID DESTRUCTION OF MAT WHILE PERFORMING ROW OPERATIONS
      TEMP = MAT

      ! PROCEDURE TO CONVERT TEMP INTO AN UPPER TRIANGULAR MATRIX WITH ALL DIAGONAL ELEMENTS EQUAL TO ONE
      DO I = 1, NDIM ! LOOP OVER ALL THE COLUMNS OF TEMP

      ! IN CASE THE ENTRY TEMP(I, I) IS ZERO, WE NEED TO FIND A GOOD PIVOT
      ! THIS PIVOT IS CHOSEN AS THE BIGGEST VALUE ON THE COLUMN J FROM TEMP(J, I) WITH J = 1, NDIM
          IF (TEMP(I, I).EQ.0.0) THEN
              BIG = I
              DO J = I, NDIM
                  IF (TEMP(J, I).GT.TEMP(BIG, I)) THEN
                      BIG = J
                  END IF
              END DO
    ! INTERCHANGE ROW I WITH BIG FOR BOTH TEMP() AND INVMAT() MATRICES
              IF (BIG.EQ.I) THEN
                  PRINT *, "SINGULAR MATRIX"
                  RETURN
              ELSE
                  DO K = 1, NDIM
                      DUMMY = TEMP(I, K)
                      TEMP(I, K) = TEMP(BIG, K)
                      TEMP(BIG, K) = DUMMY
                  
                      DUMMY = INVMAT(I, K)
                      INVMAT(I, K) = INVMAT(BIG, K)
                      INVMAT(BIG, K) = DUMMY
                  END DO
              END IF
          END IF
   
          ! DIVIDE ALL ENTRIES IN LINE I FROM TEMP(I,J) BY THE PIVOT VALUE TEMP(I,I);
          ! SAME OPERATION FOR INVMAT MATRIX
          PIVOT = TEMP(I, I)
          DO J = 1, NDIM
              TEMP(I, J) = TEMP(I, J)/PIVOT
              INVMAT(I,J) = INVMAT(I,J)/PIVOT
          END DO
          
          ! MAKE ZERO ALL ENTRIES IN THE COLUMN I BELOW THE PIVOT THROUGH ROW OPERATIONS; SIMILAR OPERATIONS FOR INVMAT
          DO J = I+1, NDIM
              FACTOR = TEMP(J, I)
              DO K = 1, NDIM
                  TEMP(J,K) = TEMP(J,K) - FACTOR*TEMP(I,K)
                  INVMAT(J,K) = INVMAT(J,K) - FACTOR*INVMAT(I,K)
              END DO
          END DO
          
      END DO

      ! MAKE THE UPPER DIAGONAL ELEMENTS OF TEMP ZERO THROUGH ROW OPERATIONS; SIMILAR OPERATIONS FOR INVMAT
      DO I = 1,NDIM - 1
          DO J = I+1, NDIM
              FACTOR = TEMP(I, J)
              DO L = 1, NDIM
                  TEMP(I,L) = TEMP(I,L) - FACTOR*TEMP(J,L)
                  INVMAT(I,L) = INVMAT(I,L) - FACTOR*INVMAT(J,L)
              END DO
          END DO
      END DO
      
      RETURN
	END SUBROUTINE INVERSE
! ---------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE INNERPROD(A, B, A_DOT_B, NDIM)
      ! SUBROUTINE TO CALCULATE THE INNER OR DOT PRODUCT OF TWO VECTORS EACH CONTAINING NDIM ELEMENTS
      
      INCLUDE 'ABA_PARAM.INC'

      INTEGER NDIM

      REAL*8 A(NDIM), B(NDIM), A_DOT_B
      
      IF (SIZE(A).EQ.SIZE(B)) THEN
          
          A_DOT_B = 0.0
          DO I = 1, NDIM
              A_DOT_B = A_DOT_B + A(I)*B(I)
          END DO
          
      ELSE
          
          WRITE(6,*) 'VECTOR DIMENSIONS DO NOT MATCH - INPUT ERROR'
      
      END IF

      RETURN
      END SUBROUTINE INNERPROD
! ---------------------------------------------------------------------------------------------------------------------------------
      
      SUBROUTINE OUTERPROD(A, B, A_DYAD_B, NDIM_A, NDIM_B)
      ! SUBROUTINE TO CALCULATE THE OUTER OR DYADIC PRODUCT OF TWO VECTORS CONTAINING NDIM_A AND NDIM_B ELEMENTS RESPECTIVELY
      
      INCLUDE 'ABA_PARAM.INC'

      INTEGER NDIM_A, NDIM_B
      
      REAL*8 A(NDIM_A), B(NDIM_B), A_DYAD_B(NDIM_A, NDIM_B)
          
      DO I = 1, NDIM_A
          DO J = 1, NDIM_B
              A_DYAD_B(I, J) = A(I)*B(J)
          END DO
      END DO
      
      RETURN
      END SUBROUTINE OUTERPROD
! ---------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE UHARD(SYIELD,HARD,EQPLAS,EQPLASRT,TIME,DTIME,TEMP,
     1     DTEMP,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,CMNAME,NSTATV,
     2     STATEV,NUMFIELDV,PREDEF,DPRED,NUMPROPS,PROPS)

      INCLUDE 'ABA_PARAM.INC'

      REAL*8 HARD(3), STATEV(NSTATV), TIME(*),
     1 PREDEF(NUMFIELDV), DPRED(1), PROPS(NUMPROPS)

      REAL*8 K, N_EXP, EQPLAS, EPS0

!     SPECIFY HARDENING PROPERTIES

      K = PROPS(3)
      EPS0 = PROPS(4)
      N_EXP = PROPS(5)
	  
!     CURRENT YIELD STRESS AND HARD VALUE

      SYIELD = K*((EPS0 + EQPLAS)**N_EXP)
      HARD(1) = K*N_EXP*((EPS0 + EQPLAS)**(N_EXP - 1.0))

      RETURN
      END SUBROUTINE UHARD
! ---------------------------------------------------------------------------------------------------------------------------------
