program main
    use iso_fortran_env, only: dp => real64
    use lattice_mod
    implicit none

    type(Lattice)  :: L
    type(UnitCell) :: U
    integer        :: ccdsto(2)
    logical        :: in_lattice
    
    U = new_unitcell(2) ! 2 dimensional unitcell
    call U%set_lvec(1, [1.0_dp, 0.0_dp]) ! Lattice vector 1 = [1, 0]
    call U%set_lvec(2, [0.0_dp, 1.0_dp]) ! Lattice vector 2 = [0, 1]
    call U%make_rvecs() ! Make reciprocal lattice vectors

    call U%add_orbital([0.0_dp, 0.0_dp]) ! Orbital 1 at bottom left of unit cell
    call U%add_orbital([0.5_dp, 0.0_dp]) ! Orbtial 2 halfway through first lattice vector of unit cell

    L = new_lattice(U, [4, 6], [.true., .true.]) ! Make a 4 x 6 periodic in both directions lattice out of U

    print *, L%ccds_sorb_to_sind([0, 0], 1) ! Print the site index of orbital 1 in the unit cell at [0, 0]
    print *, L%ccds_sorb_to_sind([4, 2], 2) ! Print the site index of orbital 2 in the unit cell at [4, 2]

    ! Displace [4, 5] from the cell at [1, 1] to end up at ccdsto, accounting for
    ! periodicity or lack of periodicity (in_lattice is .true. if this displacement
    ! ends up in the lattice, false if not).
    call L%ccds_cdsp_to_ccds([1, 1], [4, 5], ccdsto, in_lattice)
    print *, ccdsto
    print *, in_lattice

endprogram main