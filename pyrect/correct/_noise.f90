!  ----------------------------------------------------------------------------
!  echo.correct._noise
!  ===================
!
!  Fortran routines for characterizing noise in radar data.
!
!  ----------------------------------------------------------------------------


subroutine hildebrand(power, fill_value, nr, ng, P, Q, R2)
!  ----------------------------------------------------------------------------
!  Estimate the noise floor from radar power measurements (e.g., reflectivity)
!  using the method of Hildebrand and Sekhon (1974).
!
!  Parameters
!  ----------
!  power : array, dim(nr,ng), float64
!     Radar power measurements in linear units. The array must be sorted in
!     descending power values along its first dimension.
!  fill_value : float64
!     Value indicating missing or bad data in radar power measurements.
!
!  Returns
!  -------
!
!  ----------------------------------------------------------------------------

   implicit none

   integer(kind=4), intent(in)                 :: nr, ng
   real(kind=8), intent(in), dimension(nr,ng)  :: power
   real(kind=8), intent(in)                    :: fill_value
   real(kind=8), intent(out), dimension(ng)    :: P, Q, R2

   ! Define local variables
   real(kind=8), parameter   :: atol=1.d-5
   integer(kind=4)           :: ray, gate, n
   logical, dimension(nr,ng) :: valid_gate
   real(kind=8)              :: mu, var

   ! F2PY directives
   ! f2py integer(kind=4), optional, intent(in) :: nr, ng
   ! f2py real(kind=8), intent(in)              :: power
   ! f2py real(kind=8), intent(in)              :: fill_value
   ! f2py real(kind=8), intent(out)             :: P, Q, R2
   ! f2py integer(kind=4), intent(out)          :: n

   valid_gate = abs(power - fill_value) > atol

   P = fill_value
   Q = fill_value
   R2 = fill_value

   do gate = 1, ng
      do ray = 1, nr

         ! Count sample size at fixed range
         n = count(valid_gate(ray:nr,gate))

         if (n > 1) then

            mu = sum(power(ray:nr,gate), valid_gate(ray:nr,gate)) / n
            var = sum((power(ray:nr,gate) - mu)**2, valid_gate(ray:nr,gate)) / (n - 1)

            if (mu**2 / var > 1.d0) then
               P(gate) = mu
               Q(gate) = var
               R2(gate) = power(ray,gate)
               exit
            endif

         endif

      enddo
   enddo

   return

end subroutine hildebrand
