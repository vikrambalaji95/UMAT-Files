      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     1 rpl,ddsddt,drplde,drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)

      include 'aba_param.inc'

      character*80 cmname
	  
      real*8 stress(ntens), statev(nstatv), ddsdde(ntens, ntens),
     1 ddsddt(ntens), drplde(ntens), stran(ntens), dstran(ntens),
     2 time(2), predef(1), dpred(1), props(nprops), coords(3), 
     3 drot(3, 3), dfgrd0(3, 3), dfgrd1(3, 3), jstep(4)
      
! ----------------------------------------------------------------
!     umat for isotropic von mises yield criterion (j2 plasticity)
!     linear hardening law (y = y0 + h*eqplas)
!     first order forward euler explicit algorithm
! ----------------------------------------------------------------

      real*8 eqplas, e, nu, y0, h, lambda, g, elstiff(ntens, ntens), 
     1 shydro, smises, y, yieldfn, deqpl, flow(ntens), 
     2 term1, term2, term3(ntens, ntens)

!     material properties input

      e = props(1) ! young's modulus
      nu = props(2) ! poisson's ratio
      y0 = props(3) ! initial yield stress
      h = props(4) ! linear hardening rate
      
!     recover solution dependent state variables from previous step time increment
      eqplas = statev(1) ! equivalent plastic strain

!     set up elastic stiffness matrix 

      elstiff = 0.0 ! initialise the elastic stiffness matrix to zero

      lambda = nu*e/((1.0 + nu)*(1.0 - 2.0*nu)) ! 1st lamé parameter
      g = e/(2.0*(1.0 + nu)) ! 2nd lamé parameter

      if (ndi.eq.3) then ! number of direct stress components is 3
!     the problem is either 3d or plane strain
      do i = 1, ndi
          do j = 1, ndi
              elstiff(j, i) = lambda
          end do
          elstiff(i, i) = lambda + 2.0*g
      end do
      do i = ndi+1, ntens
          elstiff(i, i) = g
      end do

      else ! number of direct stress components is 2
!     the problem is plane stress
      do i = 1, ndi 
          do j = 1, ndi
              elstiff(j, i) = nu*e/(1.0 - (nu)**2.0)
          end do
          elstiff(i, i) = e/(1.0 - (nu)**2.0)
      end do
      do i = ndi+1, ntens
          elstiff(i, i) = g
      end do

      end if

      call sinv(stress, shydro, smises, ndi, nshr) !utility routine sinv to calculate the stress invariants
      ! 1st invariant is the trace of the stress matrix which gives hydrostatic stress
      ! 2nd invariant is the von mises equivalent stress - sqrt(3*j2)

      y = y0 + h*eqplas ! linear hardening law

      yieldfn = smises - y ! yield function

      if (yieldfn.lt.0.0) then
      ! the state of stress is elastic, hence zero plastic deformation in this increment
      deqpl = 0.0

      stress = stress + matmul(elstiff, dstran)
      
      ddsdde = elstiff ! for elastic case, tangent stiffness matrix is same as the elastic stiffness matrix

      else
      ! non zero plastic deformation in this increment
      
      ! flow vector for von mises
      do i = 1, ndi ! direct components
      flow(i) = 3.0*(stress(i) - shydro)/(2.0*smises)
      end do
      do i = ndi + 1, ntens ! shear components
      flow(i) = 3.0*stress(i)/smises
      end do

      ! calculate the increment in effective plastic strain

      call innerprod(flow, matmul(elstiff,dstran), term1, ntens)
      call innerprod(flow, matmul(elstiff,flow), term2, ntens)

      deqpl = term1/(term2 + h)

      ! calculate the updated stress
      stress = stress + matmul(elstiff, dstran - deqpl*flow)

      ! evaluate the material tangent stiffness matrix ddsdde
      call outerprod(matmul(elstiff,flow), matmul(elstiff,flow), 
     1 term3, ntens, ntens)

      ddsdde = elstiff - (1.0/(term2 + h))*term3

      end if

      ! update the equivalent plastic strain and store it in the solution dependent state variable
      statev(1) = eqplas + deqpl

      return
      end subroutine umat
! ---------------------------------------------------------------------------------------------------------------------------------

      subroutine innerprod(a, b, a_dot_b, ndim)
      ! subroutine to calculate the inner or dot product of two vectors each containing ndim elements
      
      include 'aba_param.inc'

      integer ndim

      real*8 a(ndim), b(ndim), a_dot_b
      
      if (size(a).eq.size(b)) then
          
          a_dot_b = 0.0
          do i = 1, ndim
              a_dot_b = a_dot_b + a(i)*b(i)
          end do
          
      else
          
          write(6,*) 'vector dimensions do not match - input error'
      
      end if

      return
      end subroutine innerprod
! ---------------------------------------------------------------------------------------------------------------------------------
      
      subroutine outerprod(a, b, a_dyad_b, ndim_a, ndim_b)
      ! subroutine to calculate the outer or dyadic product of two vectors containing ndim_a and ndim_b elements respectively
      
      include 'aba_param.inc'

      integer ndim_a, ndim_b
      
      real*8 a(ndim_a), b(ndim_b), a_dyad_b(ndim_a, ndim_b)
          
      do i = 1, ndim_a
          do j = 1, ndim_b
              a_dyad_b(i, j) = a(i)*b(j)
          end do
      end do
      
      return
      end subroutine outerprod
! ---------------------------------------------------------------------------------------------------------------------------------
