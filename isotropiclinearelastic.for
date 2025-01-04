      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1 DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP,
     2 PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS,
     3 COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER,
     4 KSPT, KSTEP, KINC)

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME

      REAL*8 STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS, NTENS),
     1 DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS),
     2 TIME(2), PREDEF(1), DPRED(1), PROPS(NPROPS), COORDS(3),
     3 DROT(3, 3), DFGRD0(3, 3), DFGRD1(3, 3), JSTEP(4)

      REAL*8 E, NU, LAMBDA, G, ELSTIFF(NTENS, NTENS)
! ----------------------------------------------------------------
!     UMAT FOR ISOTROPIC LINEAR ELASTICITY
!     APPLIES TO BOTH 3D AND PLANE STRESS SCENARIOS
! ----------------------------------------------------------------

!     MATERIAL PROPERTIES INPUT
      E = PROPS(1) ! YOUNG'S MODULUS
      NU = PROPS(2) ! POISSON'S RATIO

!     ELASTIC STIFFNESS MATRIX

      ELSTIFF = 0.0 ! Initialise the elastic stiffness matrix to zero

      LAMBDA = (NU*E)/((1.0 + NU)*(1.0 - (2.0*NU))) ! 1st Lamé Parameter
      G = E/(2.0*(1.0 + NU)) ! 2nd Lamé Parameter

      IF (NDI.EQ.3) THEN ! Number of direct stress components is 3
!     The problem is either 3d or plane strain
      DO K1=1, NDI
          DO K2=1, NDI
              ELSTIFF(K2, K1) = LAMBDA
          END DO
          ELSTIFF(K1, K1) = LAMBDA + 2.0*G
      END DO
      DO K1 = NDI+1, NTENS
          ELSTIFF(K1, K1) = G
      END DO

      ELSE ! Number of direct stress components is 2
!     The problem is plane stress
      DO K1=1, NDI 
          DO K2=1, NDI
              ELSTIFF(K2, K1) = NU*E/(1.0 - (NU)**2.0)
          END DO
          ELSTIFF(K1, K1) = E/(1.0 - (NU)**2.0)
      END DO
      DO K1 = NDI+1, NTENS
          ELSTIFF(K1, K1) = G
      END DO

      END IF

!     CALCULATE THE STRESS USING GENERALISED HOOKE'S LAW

      STRESS = STRESS + MATMUL(ELSTIFF,DSTRAN)

!     TANGENT STIFFNESS MATRIX

      DDSDDE = ELSTIFF ! For linear elastic case, tangent stiffness matrix is same as the elastic stiffness matrix

      RETURN
      END
