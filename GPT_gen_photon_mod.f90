module gen_photon_mod
use cons
use random
use voigt_mod
use mpi
use interpolate_mod
implicit none
public

public gen_photon
public gen_photon_flat
public gen_photon_AGN_con
public gen_photon_delta
public gen_photon_Gaussian
public initial_photon

contains

subroutine initial_photon(photon)
! isotropic emission
type(photon_type) :: photon
real(kind=rkd) cost,phi
real(kind=rkd) :: cosp, sinp, sint

	! initial wavevector
	cost = 2.d0*rand_number() - 1.d0
	phi = 2.d0*pi*rand_number()
	sint = sqrt(1.d0 - cost**2)
	cosp = cos(phi)
	sinp = sin(phi)
	photon%kx = sint*cosp
	photon%ky = sint*sinp
	photon%kz = cost
	photon%esc = .false.

	! For Polarization
	photon%mx = cost*cosp
	photon%my = cost*sinp
	photon%mz = -sint
	photon%nx = -sinp
	photon%ny = cosp
	photon%nz = 0.d0

	! Initial Stokes Parameter
	photon%I = 1.d0
	photon%Q = 0.d0
	photon%U = 0.d0
	photon%V = 0.d0

	photon%E3 = 1.0d0
	photon%vel1 = 0.d0
	photon%vel2 = 0.d0

	photon%weight = 1.d0
	photon%NS = 0
	photon%NS_K = 0
	photon%NS_H = 0
	photon%NS_D = 0
	photon%path = 0.d0
	photon%nclump = 0
	photon%clump = .false.
	photon%tau_atom = 0.d0
	photon%tau_dust = 0.d0

end subroutine initial_photon

subroutine gen_photon_Gaussian(photon, v_emit)
use interpolate_mod
use random
implicit none
type(photon_type) :: photon
real(kind=rkd), intent(in) :: v_emit
real(kind=rkd) vx,vy,vz
real(kind=rkd) vel
double precision, allocatable :: radius(:), emiss(:)
double precision :: radius_max, r0, px, p_max
integer :: ii

! Read emissivity data
integer :: n_points
	call read_data_file('/home/jin/CIV_Photoionization/CIV_emissivity.txt', radius, emiss, n_points)
	radius_max = 100.0d0

	! Normalize emissivity
	total_emiss = sum(emiss)
	emiss = emiss / total_emiss
	p_max = maxval(emiss)

	! Select random radius using rejection sampling
	do
	    r1,r2 = rand_number(2)
	    r0 = radius_max *r1
	    px = find_y(r0, radius, emiss, n_points)
	    if (r2 <= px / p_max) exit
	end do

	! Convert radius to coordinates
	double precision :: theta, phi, x_i, y_i, z_i
	theta = acos(2.d0 * rand_number() - 1.d0)
	phi = 2.d0 * pi * rand_number()
	x_i = r0 * sin(theta) * cos(phi) *kpc
	y_i = r0 * sin(theta) * sin(phi) *kpc
	z_i = r0 * cos(theta) * kpc

	! Assign position
	photon%x = x_i
	photon%y = y_i
	photon%z = z_i
	photon%ix = grid%N_X / 2 + 1
	photon%iy = grid%N_Y / 2 + 1
	photon%iz = grid%N_Z / 2 + 1

	photon%x_s = photon%x
	photon%y_s = photon%y
	photon%z_s = photon%z

	! Random velocity in emission region
	vel = rand_gauss() * v_emit

	! Initial wavelength (frequency)
	real(kind=rkd) :: temp
	temp = rand_number()
	if (temp <= atom%f12_K / atom%f12) then
	    photon%nu = atom%nuK / (1.d0 - vel / c)
	    photon%line = 1
	else
	    photon%nu = atom%nuH / (1.d0 - vel / c)
	    photon%line = 2
	end if

call initial_photon(photon)

end subroutine gen_photon_Gaussian



end module gen_photon_mod

