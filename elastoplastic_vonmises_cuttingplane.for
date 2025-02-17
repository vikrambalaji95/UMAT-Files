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
!     UMAT FOR ISOTROPIC VON MISES YIELD CRITERION (J2 PLASTICITY)
!     SWIFT HARDENING LAW - Y = K*(EPS0 + EQPLAS)^N
!     SEMI-EXPLICIT CUTTING PLANE PROJECTION ALGORITHM
! ----------------------------------------------------------------
      
      REAL*8 E, NU, K, EPS0, N, EQPLAS, LAMBDA, G, SMISES, SHYDRO,
     1 ELSTIFF(NTENS, NTENS), YIELDFN, HARD(3), DGAMMA, DDGAMMA,
     2 FLOW(NTENS), TERM1, TERM2(NTENS, NTENS)
      
      PARAMETER(TOLER = 1.D-6)

!     MATERIAL PROPERTIES INPUT

      E = PROPS(1)
      NU = PROPS(2)
      K = PROPS(3)
      EPS0 = PROPS(4)
      N = PROPS(5)
    
!     RECOVER SOLUTION DEPENDENT STATE VARIABLES FROM PREVIOUS STEP TIME INCREMENT
      EQPLAS = STATEV(1) ! Equivalent plastic strain

!     SET UP ELASTIC STIFFNESS MATRIX
      
      LAMBDA = NU*E/((1.0 + NU)*(1.0 - 2.0*NU)) ! 1st Lamé Parameter
      G = E/(2.0*(1.0 + NU)) ! 2nd Lamé Parameter

      ELSTIFF = 0.0 ! Initialise the elastic stiffness matrix to zero

      IF (NDI.EQ.3) THEN ! Number of direct stress components is 3
!     The problem is either 3d or plane strain
      DO I = 1, NDI
          DO J = 1, NDI
              ELSTIFF(J, I) = LAMBDA
          END DO
          ELSTIFF(I, I) = LAMBDA + 2.0*G
      END DO
      DO I = NDI+1, NTENS
          ELSTIFF(I, I) = G
      END DO

      ELSE ! Number of direct stress components is 2
!     The problem is plane stress
      DO I = 1, NDI 
          DO J = 1, NDI
              ELSTIFF(J, I) = NU*E/(1.0 - (NU)**2.0)
          END DO
          ELSTIFF(I, I) = E/(1.0 - (NU)**2.0)
      END DO
      DO I = NDI+1, NTENS
          ELSTIFF(I, I) = G
      END DO

      END IF

!     CALCULATE PREDICTOR/TRIAL STRESS ASSUMING ELASTIC BEHAVIOUR

      STRESS = STRESS + MATMUL(ELSTIFF, DSTRAN)
      
      CALL SINV(STRESS, SHYDRO, SMISES, NDI, NSHR) !Utility routine SINV to calculate the stress invariants
      ! 1st invariant is the trace of the stress matrix which gives hydrostatic stress
      ! 2nd invariant is the von mises equivalent stress - sqrt(3*J2)

      CALL UHARD(SYIELD,HARD,EQPLAS,EQPLASRT,TIME,DTIME,TEMP,
     1     DTEMP,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,CMNAME,NSTATV,
     2     STATEV,NUMFIELDV,PREDEF,DPRED,NUMPROPS,PROPS)

      YIELDFN = SMISES - SYIELD
      
      IF (YIELDFN.LT.0.0) THEN
       
!     The state of stress is elastic
!     Predicted value of stress is equal to the true stress
!     Set elastic stiffness matrix as material tangent stiffness matrix
      DDSDDE = ELSTIFF

      ELSE

!     Plastic deformation has commenced
!     Solve for equivalent plastic strain increment using the cutting plane procedure
          
      DGAMMA = 0.0 ! Initialise the increment in consistency parameter to zero
      
      DO WHILE (ABS(YIELDFN).GE.TOLER) ! The iterations last until the value of yield function becomes zero
          
      ! Flow vector for von mises
      DO I = 1, NDI
      FLOW(I) = 3.0*(STRESS(I) - SHYDRO)/(2.0*SMISES)
      END DO
      DO I = NDI + 1, NTENS
      FLOW(I) = 3.0*STRESS(I)/SMISES
      END DO
      
      CALL INNERPROD(FLOW, MATMUL(ELSTIFF, FLOW), TERM1, NTENS)
      
      DDGAMMA = YIELDFN/(TERM1 + HARD(1)) ! Current iterative increment in effective plastic strain 
      
      DGAMMA = DGAMMA + DDGAMMA ! Append increment in effective plastic strain
      
      STRESS = STRESS - MATMUL(ELSTIFF, DDGAMMA*FLOW) 
      
      CALL SINV(STRESS, SHYDRO, SMISES, NDI, NSHR)
      
      CALL UHARD(SYIELD,HARD,EQPLAS + DGAMMA,EQPLASRT,TIME,DTIME,TEMP,
     1     DTEMP,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,CMNAME,NSTATV,
     2     STATEV,NUMFIELDV,PREDEF,DPRED,NUMPROPS,PROPS)
      
      YIELDFN = SMISES - SYIELD
      
      END DO
      
      ! Evaluate the material tangent stiffness matrix DDSDDE
      CALL OUTERPROD(MATMUL(ELSTIFF,FLOW), MATMUL(ELSTIFF,FLOW), 
     1 TERM2, NTENS, NTENS)

      DDSDDE = ELSTIFF - (1.0/(TERM1 + HARD(1)))*TERM2
      
      END IF
      
      ! Update the equivalent plastic strain and store it in the solution dependent state variable
      STATEV(1) = EQPLAS + DGAMMA

      RETURN
      END SUBROUTINE UMAT     
! ---------------------------------------------------------------------------------------------------------------------------------

      SUBROUTINE INNERPROD(A, B, A_DOT_B, NDIM)
      ! Subroutine to calculate the inner or dot product of two vectors each containing NDIM elements
      
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
      ! Subroutine to calculate the outer or dyadic product of two vectors containing NDIM_A and NDIM_B elements respectively
      
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

      REAL*8 K, N, EQPLAS, EPS0

!     SPECIFY HARDENING PROPERTIES

      K = PROPS(3)
      EPS0 = PROPS(4)
      N = PROPS(5)
	  
!     CURRENT YIELD STRESS AND HARD VALUE

      SYIELD = K*((EPS0 + EQPLAS)**N)
      HARD(1) = K*N*((EPS0 + EQPLAS)**(N - 1.0))

      RETURN
      END SUBROUTINE UHARD
! ---------------------------------------------------------------------------------------------------------------------------------