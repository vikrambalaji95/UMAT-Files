      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
	  
      REAL*8 STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS, NTENS),
     1 DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS),
     2 TIME(2), PREDEF(1), DPRED(1), PROPS(NPROPS), COORDS(3), 
     3 DROT(3, 3), DFGRD0(3, 3), DFGRD1(3, 3), JSTEP(4)
      
! ----------------------------------------------------------------
!     UMAT FOR ISOTROPIC VON MISES YIELD CRITERION (J2 PLASTICITY)
!     LINEAR HARDENING LAW (Y = Y0 + H*EQPLAS)
!     FIRST ORDER FORWARD EULER EXPLICIT ALGORITHM
! ----------------------------------------------------------------

      REAL*8 EQPLAS, E, NU, Y0, H, LAMBDA, G, ELSTIFF(NTENS, NTENS), 
     1 SHYDRO, SMISES, Y, YIELDFN, DEQPL, FLOW(NTENS), 
     2 TERM1, TERM2, TERM3(NTENS, NTENS)

!     MATERIAL PROPERTIES INPUT

      E = PROPS(1) ! Young's Modulus
      NU = PROPS(2) ! Poisson's Ratio
      Y0 = PROPS(3) ! Initial yield stress
      H = PROPS(4) ! Linear hardening rate
      
!     RECOVER SOLUTION DEPENDENT STATE VARIABLES FROM PREVIOUS STEP TIME INCREMENT
      EQPLAS = STATEV(1) ! Equivalent plastic strain

      LAMBDA = NU*E/((1.0 + NU)*(1.0 - 2.0*NU)) ! 1st Lamé Parameter
      G = E/(2.0*(1.0 + NU)) ! 2nd Lamé Parameter

!     SET UP ELASTIC STIFFNESS MATRIX 

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

      CALL SINV(STRESS, SHYDRO, SMISES, NDI, NSHR) !Utility routine SINV to calculate the stress invariants
      ! 1st invariant is the trace of the stress matrix which gives hydrostatic stress
      ! 2nd invariant is the von mises equivalent stress - sqrt(3*J2)

      Y = Y0 + H*EQPLAS ! Linear hardening law

      YIELDFN = SMISES - Y ! Yield function

      IF (YIELDFN.LT.0.0) THEN
      ! The state of stress is elastic, hence zero plastic deformation in this increment
      DEQPL = 0.0

      STRESS = STRESS + MATMUL(ELSTIFF, DSTRAN)
      
      DDSDDE = ELSTIFF ! For elastic case, tangent stiffness matrix is same as the elastic stiffness matrix

      ELSE
      ! Non zero plastic deformation in this increment
      
      ! Flow vector for von mises
      DO I = 1, NDI ! Direct components
      FLOW(I) = 3.0*(STRESS(I) - SHYDRO)/(2.0*SMISES)
      END DO
      DO I = NDI + 1, NTENS ! Shear components
      FLOW(I) = 3.0*STRESS(I)/SMISES
      END DO

      ! Calculate the increment in effective plastic strain

      CALL INNERPROD(FLOW, MATMUL(ELSTIFF,DSTRAN), TERM1, NTENS)
      CALL INNERPROD(FLOW, MATMUL(ELSTIFF,FLOW), TERM2, NTENS)

      DEQPL = TERM1/(TERM2 + H)

      ! Calculate the updated stress
      STRESS = STRESS + MATMUL(ELSTIFF, DSTRAN - DEQPL*FLOW)

      ! Evaluate the material tangent stiffness matrix DDSDDE
      CALL OUTERPROD(MATMUL(ELSTIFF,FLOW), MATMUL(ELSTIFF,FLOW), 
     1 TERM3, NTENS, NTENS)

      DDSDDE = ELSTIFF - (1.0/(TERM2 + H))*TERM3

      END IF

      ! Update the equivalent plastic strain and store it in the solution dependent state variable
      STATEV(1) = EQPLAS + DEQPL

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