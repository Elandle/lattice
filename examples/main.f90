program main
    use iso_fortran_env, only: dp => real64
    use lattice_mod
    implicit none


    type(Lattice) :: Lat
    type(UnitCell) :: U

    integer :: r(2)
    logical :: in_lattice


    U = new_unitcell(2)
    call U%set_lattice_vector(1, [1.0_dp, 0.0_dp])
    call U%set_lattice_vector(2, [0.0_dp, 1.0_dp])
    call U%add_orbital([0.0_dp, 0.0_dp])
    call U%add_orbital([0.5_dp, 0.0_dp])


    Lat = new_lattice(U, [3, 3], [.false., .false.])

    print *, Lat%cellindx_from_siteindx(9)
    print *, Lat%cellindx_from_coords([1, 2])

    call Lat%cell_displacement([1, 0], [3, 2], r, in_lattice)
    print *, r
    print *, in_lattice



endprogram main