! Module: sweeps.f90


subroutine mean_fill(input, sweep_start, sweep_end, ray_window, gate_window, &
                     min_sample, fill_value, ns, nr, ng)

   implicit none

   integer(kind=4), intent(in)                   :: ns, nr, ng
   integer(kind=4), intent(in), dimension(ns)    :: sweep_start, sweep_end
   integer(kind=4), intent(in)                   :: ray_window, gate_window
   integer(kind=4), intent(in)                   :: min_sample
   real(kind=8), intent(in)                      :: fill_value
   real(kind=8), intent(inout), dimension(nr,ng) :: input



!  Define local variables =====================================================

   logical, dimension(nr,ng) :: is_valid
   integer(kind=4)            :: N_tmp
   integer(kind=4)            :: r, g, s, r0, rn, g0, gn, s0, sn

!  ============================================================================


!  F2PY directives ============================================================

!  ============================================================================

!  Create boolean arrays
   is_valid = input /= fill_value

!  Loop over all gates
   do g = 1, ng

!     Parse stencil of gates within window
      g0 = max(1, g - gate_window / 2)
      gn = min(ng, g + gate_window / 2)

!     Loop over all sweeps
      do s = 1, ns

!        Parse sweep start and end indices
         s0 = sweep_start(s)
         sn = sweep_end(s)

!        Loop over all rays in sweep
         do r = s0, sn

!        Only need to fill in gates that are missing
         if (.not. is_valid(r,g)) then

!        Parse stencil of rays within window
         r0 = max(1, r - ray_window / 2)
         rn = min(sn, r + ray_window / 2)

!        Compute sample size within the 2-D window
         N_tmp = count(is_valid(r0:rn,g0:gn))

!        Only compute the mean within the 2-D window if the sample size is
!        greater or equal to the minimum sample size specified
         if (N_tmp >= min_sample) then
         input(r,g) = sum(input(r0:rn,g0:gn), is_valid(r0:rn,g0:gn)) / N_tmp
         endif

         endif

         enddo
      enddo
   enddo

   return

end subroutine mean_fill


subroutine remove_salt(input, sweep_start, sweep_end, ray_window, &
                       gate_window, min_sample, rays_wrap, fill_value, &
                       ns, nr, ng)

   implicit none

   integer(kind=4), intent(in)                   :: ns, nr, ng
   integer(kind=4), intent(in), dimension(ns)    :: sweep_start, sweep_end
   integer(kind=4), intent(in)                   :: ray_window, gate_window
   integer(kind=4), intent(in)                   :: min_sample
   logical, intent(in)                           :: rays_wrap
   real(kind=8), intent(in)                      :: fill_value
   real(kind=8), intent(inout), dimension(nr,ng) :: input



!  Define local variables =====================================================

   logical, dimension(nr,ng) :: is_valid
   integer(kind=4)           :: sample
   integer(kind=4)           :: r, g, s, r0, rn, g0, gn, s0, sn

!  ============================================================================


!  F2PY directives ============================================================

!  ============================================================================

!  Create boolean arrays
   is_valid = input /= fill_value

!  TODO: Address rays wrap condition (e.g., PPI scans)

!  Loop over all gates
   do g = 1, ng

!     Parse stencil of gates within window
      g0 = max(1, g - gate_window / 2)
      gn = min(ng, g + gate_window / 2)

!     Loop over all sweeps
      do s = 1, ns

!        Parse sweep start and end indices
         s0 = sweep_start(s)
         sn = sweep_end(s)

!        Loop over all rays in sweep
         do r = s0, sn

!        Only need to consider gates (pixels) that are valid, i.e. not missing
         if (is_valid(r,g)) then

!        Parse stencil of rays within window
         r0 = max(1, r - ray_window / 2)
         rn = min(sn, r + ray_window / 2)

!        Compute sample size within the 2-D window
         sample = count(is_valid(r0:rn,g0:gn))

!        If the number of valid points surrounding the gate (pixel) is less
!        than the minimum sample specified, it must be noise
         if (sample < min_sample) then
            input(r,g) = fill_value
         endif

         endif

         enddo
      enddo
   enddo

   return

end subroutine remove_salt
