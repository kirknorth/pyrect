! ----------------------------------------------------------------------------
! echo.util._texture
! ==================
!
! Fortran routines for computing texture fields. A texture field is defined as
! the standard deviation of a radar measurement within a 1-D or 2-D window
! centered around a radar gate.
!
! ----------------------------------------------------------------------------

! Necessary and/or potential future improvements to _texture submodule:
!
! * Properly handle contiguous sweep volumes, e.g., 360 deg PPI volumes.


subroutine compute_texture(field, sweep_start, sweep_end, ray_window, &
                          gate_window, rays_wrap_around, fill_value, &
                          debug, verbose, nr, ng, ns, sigma, sample_size)
! ----------------------------------------------------------------------------
! Compute the texture field (standard deviation) of a radar measurement field
! within a specified 1-D or 2-D texture window.
!
! Parameters
! ----------
! field : array, dim(nr,ng), float64
!   Radar data field.
! sweep_start : array, dim(ns), int32
!   Sweep start ray indices indicating the start of a radar sweep.
! sweep_end : array, dim(ns), int32
!   Sweep end ray indices indicating the end of a radar sweep.
! ray_window : int32
!   Number of rays in the texture window.
! gate_window : int32
!   Number of range gates in the texture window.
! rays_wrap_around : logical
!   True to indicate that all sweeps have contiguous rays, e.g., 360 deg PPI
!   radar volumes.
! fill_value : float64
!   Value indicating missing or bad data in field.
! debug : logical
!   True to print debugging information, False to suppress.
! verbose : logical
!   True to print relevant information, False to suppress.
!
! Returns
! -------
! sigma : array, dim(nr,ng), float64
!   Radar data texture field.
! sample_size : array, dim(nr,ng), int32
!   Number of radar gates available in the texture window to compute the
!   texture field.
!
! ----------------------------------------------------------------------------

  implicit none

  integer(kind=4), intent(in)                    :: nr, ng, ns
  real(kind=8), intent(in), dimension(nr,ng)     :: field
  integer(kind=4), intent(in), dimension(ns)     :: sweep_start, sweep_end
  integer(kind=4), intent(in)                    :: ray_window, gate_window
  logical, intent(in)                            :: rays_wrap_around
  real(kind=8), intent(in)                       :: fill_value
  logical, intent(in)                            :: debug, verbose
  real(kind=8), intent(out), dimension(nr,ng)    :: sigma
  integer(kind=4), intent(out), dimension(nr,ng) :: sample_size

! Define local variables
  real(kind=8), parameter    :: atol=1.d-5
  integer(kind=4)            :: sweep, gate, ray, ray_start, ray_stop
  integer(kind=4)            :: r0, rf, g0, gf
  integer(kind=4)            :: n
  real(kind=8)               :: mu, Var
  logical, dimension(nr,ng)  :: valid_gate

! F2PY directives
! f2py integer(kind=4), optional, intent(in) :: nr, ng, ns
! f2py real(kind=8), intent(in)              :: field
! f2py integer(kind=4), intent(in)           :: sweep_start, sweep_end
! f2py integer(kind=4), intent(in)           :: ray_window, gate_window
! f2py logical, intent(in)                   :: rays_wrap_around
! f2py real(kind=8), intent(in)              :: fill_value
! f2py logical, intent(in)                   :: debug, verbose
! f2py real(kind=8), intent(out)             :: sigma
! f2py integer(kind=4), intent(out)          :: sample_size

  sigma = fill_value
  valid_gate = abs(field - fill_value) > atol

  do gate = 1, ng

    ! Gates along a ray are by definition not contiguous so no wrapping
    ! conditions have to be checked
    g0 = max(1, gate - gate_window / 2)
    gf = min(ng, gate + gate_window / 2)

    do sweep = 1, ns

    ray_start = sweep_start(sweep) + 1
    ray_stop = sweep_end(sweep) + 1

    do ray = ray_start, ray_stop

      r0 = ray - ray_window / 2
      rf = ray + ray_window / 2

      ! Check for condition where current ray is close to a sweep boundary
      if (r0 < ray_start .or. rf > ray_stop) then

        ! If sweep rays are contiguous then special care is needed
        if (rays_wrap_around) then
          stop
        else
         if (r0 < ray_start) then
            r0 = ray_start
          endif
          if (rf > ray_stop) then
            rf = ray_stop
          endif
        endif

      endif

      ! Compute sample size within texture window
      n = count(valid_gate(r0:rf,g0:gf))
      sample_size(ray,gate) = n

      if (n > 1 .and. valid_gate(ray,gate)) then

      mu = sum(field(r0:rf,g0:gf), valid_gate(r0:rf,g0:gf)) / n
      Var = sum((field(r0:rf,g0:gf) - mu)**2, valid_gate(r0:rf,g0:gf)) / (n - 1)
      sigma(ray,gate) = sqrt(Var)

      endif

    enddo
    enddo
  enddo

  return

end subroutine compute_texture
