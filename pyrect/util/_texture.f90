! -----------------------------------------------------------------------------
! echo.util._texture
! ==================
!
! Fortran routines for computing textures. A texture is defined as the standard
! deviation of a radar measurement within a 1-D or 2-D window centered around a
! radar gate.
!
! -----------------------------------------------------------------------------

! Necessary and/or potential future improvements to _texture submodule:
!
! * Properly handle contiguous sweep volumes, e.g., 360 deg PPI volumes.


subroutine compute_texture(data, sweep_start, sweep_end, ray_size, &
                           gate_size, rays_wrap_around, fill_value, &
                           debug, verbose, nr, ng, ns, std, sample_size)
! -----------------------------------------------------------------------------
! Compute texture of a radar measurement.
!
! Parameters
! ----------
! data : array, dim(nr, ng), float64
!   Radar data.
! sweep_start : array, dim(ns), int32
!   Sweep start ray indices indicating the start of a radar sweep.
! sweep_end : array, dim(ns), int32
!   Sweep end ray indices indicating the end of a radar sweep.
! ray_size : int32
!   Number of rays in the texture filter.
! gate_size : int32
!   Number of range gates in the texture filter.
! rays_wrap_around : logical
!   True if rays are contiguous in all sweeps.
! fill_value : float64
!   Value indicating missing or bad data in field.
! debug : logical
!   True to print debugging information, False to suppress.
! verbose : logical
!   True to print relevant information, False to suppress.
!
! Returns
! -------
! std : array, dim(nr, ng), float64
!   Texture of radar measurement.
! sample_size : array, dim(nr, ng), int32
!   Sample size in texture filter.
!
! -----------------------------------------------------------------------------
  implicit none

  integer(kind=4), intent(in)                     :: nr, ng, ns
  real(kind=8), intent(in), dimension(nr, ng)     :: data
  integer(kind=4), intent(in), dimension(ns)      :: sweep_start, sweep_end
  integer(kind=4), intent(in)                     :: ray_size, gate_size
  logical, intent(in)                             :: rays_wrap_around
  real(kind=8), intent(in)                        :: fill_value
  logical, intent(in)                             :: debug, verbose
  real(kind=8), intent(out), dimension(nr, ng)    :: std
  integer(kind=4), intent(out), dimension(nr, ng) :: sample_size

! Define local variables
  real(kind=8), parameter    :: atol=1.d-5
  integer(kind=4)            :: sweep, gate, ray, ray_start, ray_stop
  integer(kind=4)            :: r0, rf, g0, gf
  integer(kind=4)            :: n
  real(kind=8)               :: mu, var
  logical, dimension(nr, ng) :: valid_gate

! F2PY directives
! f2py integer(kind=4), optional, intent(in) :: nr, ng, ns
! f2py real(kind=8), intent(in)              :: data
! f2py integer(kind=4), intent(in)           :: sweep_start, sweep_end
! f2py integer(kind=4), intent(in)           :: ray_size, gate_size
! f2py logical, intent(in)                   :: rays_wrap_around
! f2py real(kind=8), intent(in)              :: fill_value
! f2py logical, intent(in)                   :: debug, verbose
! f2py real(kind=8), intent(out)             :: std
! f2py integer(kind=4), intent(out)          :: sample_size

  std = fill_value
  valid_gate = abs(data - fill_value) > atol

  do gate = 1, ng

    ! Gates along a ray are by definition not contiguous so no wrapping
    ! conditions have to be checked
    g0 = max(1, gate - gate_size / 2)
    gf = min(ng, gate + gate_size / 2)

    do sweep = 1, ns

    ray_start = sweep_start(sweep) + 1
    ray_stop = sweep_end(sweep) + 1

    do ray = ray_start, ray_stop

      r0 = ray - ray_size / 2
      rf = ray + ray_size / 2

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

      ! Compute sample size within texture filter
      n = count(valid_gate(r0:rf, g0:gf))
      sample_size(ray, gate) = n

      if (n > 1 .and. valid_gate(ray,gate)) then

      mu = sum(data(r0:rf, g0:gf), valid_gate(r0:rf, g0:gf)) / n
      var = sum((data(r0:rf, g0:gf) - mu)**2.d0, valid_gate(r0:rf, g0:gf)) / (n - 1)
      std(ray, gate) = sqrt(var)

      endif

    enddo
    enddo
  enddo

  return

end subroutine compute_texture
