!IRAND = MOD(IRAND * B + 1, A)
real*8 function random1() result(rand)
   integer:: ig=538247
   common /random_dat/ ig
! local params
   integer m,m1,mult
   parameter (m=100000000,m1=10000,mult=31415821)
!  how many random numbers to store locally
   integer maxloc
   parameter (maxloc = 20)
! local vars
   logical lnew
   integer irand
   save lnew,irand
!
   integer irandh,irandl,multh,multl,iused,i
   integer isave(maxloc)
   real r,rsave(maxloc)
   save iused,isave,rsave
! data
   data irand/0/
   data lnew/.true./
! begin
   if (lnew) then
      lnew = .false.
      irand = mod(iabs(ig),m)
      iused = maxloc
   endif

   if (iused .eq. maxloc) then
!     have to generate new numbers which we store
!     into isave and rsave
!
      do i=1,maxloc
!*****multiply irand by mult, but take into account that overflow must
!*****be discarded, and do not generate an error.
!
         irandh = irand/m1
         irandl = mod(irand,m1)
         multh = mult/m1
         multl = mod(mult,m1)
!
         irand = mod(irandh*multl+irandl*multh,m1)*m1 + irandl*multl
         irand = mod(irand+1,m)
!
!*****convert irand to a real random number between 0 and 1.
!
         r = real(irand/10)*10/real(m)
         if ((r .le. 0.0) .or. (r .gt. 1.0)) r = 0.0
         rsave(i) = r
         isave(i) = irand
      end do
      iused = 0
   endif

!  now give out a number
   iused = iused + 1
   rand = rsave(iused)
   ig   = isave(iused)
end function
!
!
subroutine initRandom1(ig0)
   integer:: ig0, ig
   common /random_dat/ ig
   ig = ig0
end subroutine
!
!"normal" generates a random number from a normal Gaussian
! distribution with a mean of zero and a variance of one
function normal ()
   implicit none
   real*8 random,v1,v2,rsq
   real*8 factor,store,normal
   logical compute
   save compute,store
   data compute  / .true. /
!
!  get a pair of random values from the distribution
   if (compute) then
      rsq = 0.
      do while(rsq >= 1.0d0)
         v1 = 2.0d0 * random () - 1.0d0
         v2 = 2.0d0 * random () - 1.0d0
         rsq = v1**2 + v2**2
      end do
      factor = sqrt(-2.0d0*log(rsq)/rsq)
      store = v1 * factor
      normal = v2 * factor
      compute = .false.
!
! use the second random value computed at the last call
   else
      normal = store
      compute = .true.
   end if
   return
end function
