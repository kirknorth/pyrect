! Module: texture.f90


subroutine compute_texture(input, sweep_start, sweep_end, ray_window, &
	                       gate_window, fill_value, ns, nr, ng, sample, &
	                       texture)

	implicit none

	integer(kind=4), intent(in)                    :: ns, nr, ng
	integer(kind=4), intent(in), dimension(ns)     :: sweep_start, sweep_end
	integer(kind=4), intent(in)                    :: ray_window, gate_window
	real(kind=8), intent(in)                       :: fill_value
	real(kind=8), intent(in), dimension(nr,ng)     :: input
	integer(kind=4), intent(out), dimension(nr,ng) :: sample
	real(kind=8), intent(out), dimension(nr,ng)    :: texture

!	Define local variables ===================================================

	logical, dimension(nr, ng) :: mask
	real(kind=8)               :: N_tmp, mean_tmp, var_tmp
	integer(kind=4)            :: r, g, s, r0, rn, g0, gn, s0, sn

! 	==========================================================================


!	F2PY directives ==========================================================

	!f2py integer(kind=4), optional, intent(in) :: ns, nr, ng


!	==========================================================================
	
!	Fill texture array
	texture = fill_value

!	Create boolean array where good gates are True and bad gates are False
	mask = input /= fill_value

	do g = 1, ng

!		Parse stencil of gates within window
	    g0 = max(1, g - gate_window / 2)
	    gn = min(ng, g + gate_window / 2)

	  	do s = 1, ns
	  		
!		Parse sweep start and end indices
	  	s0 = sweep_start(s)
	  	sn = sweep_end(s)

	  	do r = s0, sn

!		Parse stencil of rays within window
		r0 = max(1, r - ray_window / 2)
		rn = min(sn, r + ray_window / 2)

!		Compute sample size within the window
		N_tmp = count(mask(r0:rn,g0:gn))
		sample(r,g) = N_tmp

!		Do not compute the texture field (standard deviation) if there are
!		no good gates
		if (N_tmp > 0) then
		mean_tmp = sum(input(r0:rn,g0:gn), mask(r0:rn,g0:gn)) / N_tmp
		var_tmp = sum(input(r0:rn,g0:gn)**2, mask(r0:rn,g0:gn)) / N_tmp - mean_tmp**2
		texture(r,g) = sqrt(var_tmp)
		endif

		enddo
	  	enddo
	enddo

	return

end subroutine compute_texture
