! Module: member.f90


subroutine conditional_all(input, pdf, bins, zero, fill_value, &
                           nr, ng, nb, P_cond)

   implicit none

   integer(kind=4), intent(in)                 :: nr, ng, nb
   real(kind=8), intent(in)                    :: zero, fill_value
   real(kind=8), intent(in), dimension(nb)     :: pdf, bins
   real(kind=8), intent(in), dimension(nr,ng)  :: input
   real(kind=8), intent(out), dimension(nr,ng) :: P_cond

!  Define local variables =====================================================

   logical, dimension(nr,ng) :: is_valid
   integer(kind=4)           :: i, r, g

!  ============================================================================

!  Fill the conditional probability array
   P_cond = fill_value

!  Determine which gates are valid
   is_valid = input /= fill_value

   do g = 1, ng
      do r = 1, nr

!        Only compute the conditional probability if the gate is valid
         if (is_valid(r,g)) then

!           Determine which bin the radar gate belongs to
            i = minloc(abs(bins - input(r,g)), dim=1)

!           Compute the conditional probability of the input data given its
!           bin and corresponding probability density function
!           Must also check for the zero probability condition
            if (pdf(i) >= zero) then
               P_cond(r,g) = pdf(i)
            else
               P_cond(r,g) = zero
            endif

         endif

      enddo
   enddo

   return

end subroutine conditional_all
