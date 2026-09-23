!
! gfortran -c -ffree-line-length-none lattice_mod.f90
!
! Todo:
! Make set_lvec accept real or integer vectors.
! Automatically construct reciprocal lattice vectors once
! all lattice vectors are set (and automatically update if changed).
! think if orbital positions should be set in terms of cooridnates or real position (or both)
! get displacement vectors between sites (accounting for boundary conditions)
!
!
!
!
module lattice_mod
    use iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
    implicit none
    private

    public :: UnitCell
    public :: Bond
    public :: Lattice
    
    public :: new_unitcell
    public :: new_bond
    public :: new_lattice
    public :: new_hopping

    integer , parameter :: dp = real64
    real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)

    type :: UnitCell
        integer  :: dim       ! Dimension of the unit cell
        integer  :: norbitals ! Number of orbitals in the unit cell
        real(dp) :: V         ! Volume of the unit cell

        real(dp), allocatable :: lvecs(:, :)             ! Lattice vectors
        real(dp), allocatable :: orbital_positions(:, :) ! Orbital positions
        real(dp), allocatable :: rvecs(:, :)             ! Reciprocal lattice vectors

        contains
            procedure set_lvec
            procedure lvec
            procedure add_orbital
            procedure make_rvecs
            procedure rvec
    endtype UnitCell

    type :: Bond
        integer :: orbto   ! Orbital index the bond is to.
        integer :: orbfr   ! Orbital index the bond is from.   

        integer, allocatable :: cdsp(:) ! Unit cell displacement.
                                        ! Measured in terms of lattice vectors.
    endtype Bond

    type :: Hopping
        complex(dp) :: tto   ! Hopping from  ---  t_to   ---> to
        complex(dp) :: tfr ! Hopping from <--- t_from  ---  to

        ! Internally store hoppings as complex numbers, then check later on if they are real or complex
        logical :: iscomplex
    endtype Hopping


    !> Type for representing lattices.
    !!
    !!
    !!
    type :: Lattice
        type(UnitCell)              :: U           !< Unit cell of the lattice.
        integer       , allocatable :: L(:)        !< Integer extent of lattice in units of primitive vectors.
                                                   !! `L(i) =` number of unit cells along axis defined by the `i`th lattice vector.
        logical       , allocatable :: periodic(:) !< Indicates which primitive vector directions are periodic. If `periodic(i) = .true.` then the `i`th direction is periodic otherwise it is not.
        integer                     :: ncells      !< Number of unit cells in the lattice.
                                                   !! Equal to the product of all entries in `L` (eg, `ncells = L(1) * L(2)` for a `dim = 2` lattice).
        integer                     :: nsites      !< Number of sites. nsites = ncells * norbitals (norbitals = number of orbitals per site).
        integer                     :: dim         !< Dimension of lattice. Can only be `1`, `2`, or `3`.
        type(Bond)    , allocatable :: bonds(:)    !< Bonds of the lattice.
        integer                     :: nbonds      !< Number of bonds in the lattice.
        type(Hopping) , allocatable :: hoppings(:)
        integer                     :: nhoppings
        logical                     :: complexhoppings


    contains
        !> @name lattice procedure naming
        !! @{
        !! Naming conventions for positioning and indexing procedures:
        !!
        !! - `ind` = index
        !! - `cds` = coordinates (in lattice vector basis)
        !! - `pos` = position (in Cartesian basis)
        !! - `orb` = orbital
        !! - `dsp` = displacement (in lattice vector basis)
        !! - `sep` = separation (in Cartesian basis)
        !! - `s` prefix = site
        !! - `c` prefix = unit cell
        !! - `fr` suffix = from
        !! - `to` suffix = to
        !!
        !! Some explanation and rationale:
        !!
        !! The inclusion of many site and cell indexing and positioning
        !! routines is meant to be exhaustive for what you might need in practice (but is not totally exhaustive).
        !! In particular, effort was made to keep as many things integer-based as possible with as little
        !! real-based inputs as possible. This is to avoid having to deal with dealing with possibly rounding real numbers
        !! or doing expensive nearest-position searches. Since sites are only located at (integer-indexed) orbitals
        !! within unit cells at integer-coordinate (in the lattice vector basis) positions, then integers
        !! are really all you need. Many times you do need some Cartesian real number positions (eg, in some Monte Carlo measurements),
        !! so routines are provided for that (eg, `sind_to_spos` gets the Cartesian position of the input site index).
        !!
        !! The indexing and positioning routines generally have the following structure: `initial_to_final`.
        !! `initial` specifies how you get the result you want (eg, an index) and `final` specifies what you want (eg, a position).
        !! Some routines have a displacement (`dsp`) input. To keep things integer-based, these are restricted to cells
        !! (so all `dsp` arguments are integer arrays). These routines all have an `in_lattice` logical argument indicating whether or not the
        !! displacement lies inside of the lattice (`.true.`) or not (`.false.`) (periodicity is accounted for).
        !> @}

        procedure :: cind_to_ccds            !> Returns the coordinates of the cell with the input index.
        procedure :: cind_to_cpos            !> Returns the position of the cell with the input index.
        procedure :: cind_sorb_to_sind       !> Returns the index of the site in the provided cell at the input orbital.
        procedure :: cind_cdsp_to_cind       !> Computes the index of the cell reached from a provided cell (by index) by a displacement, keeping track of whether the displacement is in the lattice or not.
        procedure :: cind_cdsp_sorb_to_sind  !> Computes the index of the site reached from a provided cell (by index) by a displacement and specified orbital, keeping track of whether the displacement is in the lattice or not.

        procedure :: ccds_to_cind            !> Returns the index of the cell with the input coordinates.
        procedure :: ccds_sorb_to_sind       !> Returns the index of the site within the cell at the specified coordinates and specified orbital.
        procedure :: ccds_cdsp_to_cind       !> Computes the index of the cell reached
        procedure :: ccds_cdsp_to_ccds

        procedure :: sind_to_scds
        procedure :: sind_to_spos
        procedure :: sind_to_cind
        procedure :: sind_to_oind
        procedure :: sind_to_cind_oind
        procedure :: sind_cdsp_sorb_to_sind
        procedure :: sind_cdsp_to_cind
        procedure :: sind_info

        procedure :: cds_dsp_to_cds
        procedure :: cds_to_pos_dp
        procedure :: cds_to_pos_int
        generic   :: cds_to_pos => cds_to_pos_dp, cds_to_pos_int

        procedure :: pos_to_cds

        procedure :: add_bond
        procedure :: get_bond
        procedure :: set_bond
        procedure :: deallocate_bonds
        procedure :: allocate_bonds

        procedure :: add_hopping
        procedure :: get_hopping
        procedure :: set_hopping
        procedure :: deallocate_hoppings
        procedure :: allocate_hoppings

        procedure :: ccds_ccds_to_cdsp
        procedure :: cind_cind_to_cdsp
        procedure :: sind_sind_to_cdsp
        procedure :: oind_cdsp_oind_to_sep
        procedure :: sind_sind_to_sep

        procedure :: hopping_matrix_zp
    endtype Lattice

    interface append_column
        module procedure :: append_column_dp
    endinterface append_column

    interface isclose
        module procedure :: isclose_dp
        module procedure :: isclose_zp
    endinterface isclose

    interface iscomplex
        module procedure :: iscomplex_zp
    endinterface iscomplex

    contains
        subroutine ccds_cdsp_to_ccds(self, ccdsfr, cdsp, ccdsto, in_lattice)
            class(Lattice), intent(in)  :: self
            integer       , intent(in)  :: ccdsfr(self%dim)
            integer       , intent(in)  :: cdsp(self%dim)
            integer       , intent(out) :: ccdsto(self%dim)
            logical       , intent(out) :: in_lattice

            call self%cds_dsp_to_cds(ccdsfr, cdsp, ccdsto, in_lattice)
        endsubroutine ccds_cdsp_to_ccds

        function cds_to_pos_dp(self, cds) result(pos)
            class(Lattice), intent(in) :: self
            real(dp)      , intent(in) :: cds(self%dim)

            real(dp) :: pos(self%dim)

            pos = matmul(self%U%lvecs, cds)
        endfunction cds_to_pos_dp

        function cds_to_pos_int(self, cds) result(pos)
            class(Lattice), intent(in) :: self
            integer       , intent(in) :: cds(self%dim)

            real(dp) :: pos(self%dim)

            pos = self%cds_to_pos_dp(real(cds, dp))
        endfunction cds_to_pos_int

        subroutine ccds_cdsp_to_cind(self, ccdsfr, cdsp, cindto, in_lattice)
            class(Lattice), intent(in)  :: self
            integer       , intent(in)  :: ccdsfr(self%dim)
            integer       , intent(in)  :: cdsp(self%dim)
            integer       , intent(out) :: cindto
            logical       , intent(out) :: in_lattice

            integer :: ccdsto(self%dim)

            call self%cds_dsp_to_cds(ccdsfr, cdsp, ccdsto, in_lattice)

            if (in_lattice) cindto = self%ccds_to_cind(ccdsto)
        endsubroutine ccds_cdsp_to_cind

        subroutine sind_cdsp_to_cind(self, sindfr, cdsp, cindto, in_lattice)
            class(Lattice), intent(in)  :: self
            integer       , intent(in)  :: sindfr
            integer       , intent(in)  :: cdsp(self%dim)
            integer       , intent(out) :: cindto
            logical       , intent(out) :: in_lattice

            integer :: cindfr

            cindfr = self%sind_to_cind(sindfr)

            call self%cind_cdsp_to_cind(cindfr, cdsp, cindto, in_lattice)
        endsubroutine sind_cdsp_to_cind

        subroutine cind_cdsp_sorb_to_sind(self, cindfr, cdsp, sorbto, sindto, in_lattice)
            class(Lattice), intent(in)  :: self
            integer       , intent(in)  :: cindfr
            integer       , intent(in)  :: cdsp(self%dim)
            integer       , intent(in)  :: sorbto
            integer       , intent(out) :: sindto
            logical       , intent(out) :: in_lattice

            integer :: cindto

            call self%cind_cdsp_to_cind(cindfr, cdsp, cindto, in_lattice)

            if (in_lattice) sindto = self%cind_sorb_to_sind(sorbto, cindto)
        endsubroutine cind_cdsp_sorb_to_sind

        function pos_to_cds(self, pos) result(cds)
            class(Lattice), intent(in) :: self
            real(dp)      , intent(in) :: pos(self%dim)

            real(dp) :: cds(self%dim)

            integer :: i

            if (.not. allocated(self%U%rvecs)) then
                error stop "error stop in procedure pos_to_cds from module lattice_mod: reciprocal lattice vectors have not been made before calling."
            endif

            do i = 1, self%dim
                cds(i) = dot_product(pos, self%U%rvecs(:, i)) / (2.0_dp * pi)
            enddo
        endfunction pos_to_cds

        function sind_to_oind(self, sind) result(oind)
            class(Lattice), intent(in) :: self
            integer       , intent(in) :: sind

            integer :: oind

            oind = modulo(sind - 1, self%U%norbitals) + 1
        endfunction sind_to_oind

        subroutine sind_to_cind_oind(self, sind, cind, oind)
            class(Lattice), intent(in)  :: self
            integer       , intent(in)  :: sind
            integer       , intent(out) :: cind
            integer       , intent(out) :: oind

            cind = self%sind_to_cind(sind)
            oind = self%sind_to_oind(sind)
        endsubroutine sind_to_cind_oind

        subroutine sind_info(self, sind, pos, cds, oind, cind)
            class(Lattice), intent(in) :: self
            integer       , intent(in) :: sind

            real(dp), intent(out), optional :: pos(self%dim)
            real(dp), intent(out), optional :: cds(self%dim)
            integer , intent(out), optional :: oind
            integer , intent(out), optional :: cind

            integer  :: cindtemp, oindtemp
            integer  :: ccds(self%dim)
            real(dp) :: scds(self%dim)

            call self%sind_to_cind_oind(sind, cindtemp, oindtemp)

            if (present(cind)) cind = cindtemp
            if (present(oind)) oind = oindtemp

            if (present(cds) .or. present(pos)) then
                ccds = self%cind_to_ccds(cindtemp)
                scds = real(ccds, dp) + self%U%orbital_positions(:, oindtemp)

                if (present(cds)) cds = scds
                if (present(pos)) pos = self%cds_to_pos(scds)
            endif
        endsubroutine sind_info

        function sind_to_scds(self, sind) result(scds)
            class(Lattice), intent(in) :: self
            integer       , intent(in) :: sind

            real(dp) :: scds(self%dim)
            integer  :: cind, oind
            integer  :: ccds(self%dim)

            call self%sind_to_cind_oind(sind, cind, oind)
            ccds = self%cind_to_ccds(cind)
            scds = real(ccds, dp) + self%U%orbital_positions(:, oind)
        endfunction sind_to_scds

        function sind_to_spos(self, sind) result(spos)
            class(Lattice), intent(in) :: self
            integer       , intent(in) :: sind

            real(dp) :: spos(self%dim)
            real(dp) :: scds(self%dim)

            scds = self%sind_to_scds(sind)
            spos = self%cds_to_pos(scds)
        endfunction sind_to_spos

        function cind_to_cpos(self, cind) result(cpos)
            class(Lattice), intent(in) :: self
            integer       , intent(in) :: cind

            real(dp) :: cpos(self%dim)
            integer  :: ccds(self%dim)
            real(dp) :: cds(self%dim)

            ccds = self%cind_to_ccds(cind)
            cds  = real(ccds, dp)
            cpos = self%cds_to_pos(cds)
        endfunction cind_to_cpos

        !> Adds a bond to the lattice.
        !!
        !! @param[in,out] self Lattice to add a bond to.
        !! @param[in]     B    Bond to add to lattice.
        subroutine add_bond(self, B)
            class(Lattice), intent(inout) :: self
            type(Bond)    , intent(in)    :: B

            type(Bond), allocatable :: temp(:)

            ! Never added a bond before? If so, get things ready and add this new bond.
            if (.not. allocated(self%bonds)) then
                allocate(self%bonds(1))
                self%bonds(1) = B
                self%nbonds = 1
                return
            endif

            ! Add space for 1 more bond and put this new bond there.
            allocate(temp(self%nbonds+1))
            temp(1:self%nbonds)   = self%bonds
            temp(self%nbonds + 1) = B
            call move_alloc(temp, self%bonds)
            self%nbonds = self%nbonds + 1
        endsubroutine add_bond

        subroutine add_hopping(self, t)
            class(Lattice), intent(inout) :: self
            type(Hopping)    , intent(in)    :: t

            type(Hopping), allocatable :: temp(:)

            ! Never added a hopping before? If so, get things ready and add this new hopping.
            if (.not. allocated(self%hoppings)) then
                allocate(self%hoppings(1))
                self%hoppings(1) = t
                self%nhoppings = 1
                return
            endif

            ! Add space for 1 more hopping and put this new hopping there.
            allocate(temp(self%nhoppings+1))
            temp(1:self%nhoppings)   = self%hoppings
            temp(self%nhoppings + 1) = t
            call move_alloc(temp, self%hoppings)
            self%nhoppings = self%nhoppings + 1
        endsubroutine add_hopping

        function get_bond(self, i) result(B)
            class(Lattice), intent(in) :: self
            integer       , intent(in) :: i

            type(Bond) :: B

            B = self%bonds(i)
        endfunction get_bond

        function get_hopping(self, i) result(t)
            class(Lattice), intent(in) :: self
            integer       , intent(in) :: I

            type(Hopping) :: t

            t = self%hoppings(i)
        endfunction get_hopping

        subroutine set_bond(self, i, B)
            class(Lattice), intent(inout) :: self
            integer       , intent(in)    :: i
            type(Bond)    , intent(in)    :: B

            self%bonds(i) = B
        endsubroutine set_bond

        subroutine set_hopping(self, i, t)
            class(Lattice), intent(inout) :: self
            integer       , intent(in)    :: i
            type(Hopping) , intent(in)    :: t

            self%hoppings(i) = t
        endsubroutine set_hopping

        subroutine deallocate_bonds(self)
            class(Lattice), intent(inout) :: self

            if (allocated(self%bonds)) deallocate(self%bonds)
            self%nbonds = 0
        endsubroutine deallocate_bonds

        subroutine deallocate_hoppings(self)
            class(Lattice), intent(inout) :: self

            if (allocated(self%hoppings)) deallocate(self%hoppings)
            self%nhoppings = 0
        endsubroutine deallocate_hoppings

        subroutine allocate_bonds(self, nbonds)
            class(Lattice), intent(inout) :: self
            integer       , intent(in)    :: nbonds

            call self%deallocate_bonds()
            allocate(self%bonds(nbonds))
            self%nbonds = nbonds
        endsubroutine allocate_bonds

        subroutine allocate_hoppings(self, nhoppings)
            class(Lattice), intent(inout) :: self
            integer       , intent(in)    :: nhoppings

            call self%deallocate_hoppings()
            allocate(self%hoppings(nhoppings))
            self%nhoppings = nhoppings
        endsubroutine allocate_hoppings


        subroutine sind_cdsp_sorb_to_sind(self, sindfr, cdsp, sorbto, sindto, in_lattice)
            class(Lattice), intent(in)  :: self
            integer       , intent(in)  :: sindfr
            integer       , intent(in)  :: cdsp(self%dim)
            integer       , intent(in)  :: sorbto
            integer       , intent(out) :: sindto
            logical       , intent(out) :: in_lattice

            integer :: cindfr, cindto

            cindfr = self%sind_to_cind(sindfr)
            call self%cind_cdsp_to_cind(cindfr, cdsp, cindto, in_lattice)
            if (in_lattice) sindto = self%cind_sorb_to_sind(sorbto, cindto)
        endsubroutine sind_cdsp_sorb_to_sind

        subroutine cind_cdsp_to_cind(self, cindfr, cdsp, cindto, in_lattice)
            class(Lattice), intent(in)  :: self
            integer       , intent(in)  :: cindfr
            integer       , intent(in)  :: cdsp(self%dim)
            integer       , intent(out) :: cindto
            logical       , intent(out) :: in_lattice

            integer :: ccdsfr(self%dim)
            integer :: ccdsto(self%dim)

            ccdsfr = self%cind_to_ccds(cindfr)
            call self%cds_dsp_to_cds(ccdsfr, cdsp, ccdsto, in_lattice)
            if (in_lattice) cindto = self%ccds_to_cind(ccdsto)
        endsubroutine cind_cdsp_to_cind

        
        subroutine cds_dsp_to_cds(self, cdsfr, dsp, cdsto, in_lattice)
            class(Lattice), intent(in)  :: self
            integer       , intent(in)  :: cdsfr(self%dim)
            integer       , intent(in)  :: dsp(self%dim)
            integer       , intent(out) :: cdsto(self%dim)
            logical       , intent(out) :: in_lattice

            integer :: i

            associate (dim => self%dim, L => self%L, periodic => self%periodic)

                cdsto = cdsfr + dsp
                in_lattice = .true.

                do i = 1, dim
                    if (periodic(i)) then
                        cdsto(i) = modulo(cdsto(i), L(i))
                    elseif ((cdsto(i) .lt. 0) .or. (cdsto(i) .ge. L(i))) then
                        in_lattice = .false.
                        return
                    endif
                enddo
            endassociate
        endsubroutine cds_dsp_to_cds

        function cind_to_ccds(self, cind) result(ccds)
            class(Lattice), intent(in) :: self
            integer, intent(in) :: cind

            integer :: ccds(self%dim)
            integer :: i, n

            n = cind - 1

            do i = 1, self%dim
                ccds(i) = modulo(n, self%L(i))
                n = n / self%L(i)
            enddo
        endfunction cind_to_ccds

        function new_lattice(U, L, periodic) result(Lat)
            type(UnitCell), intent(in) :: U
            integer       , intent(in) :: L(:)
            logical       , intent(in) :: periodic(:)

            type(Lattice) :: Lat

            integer :: i

            ! Check that the provided L extent matches the dimension of the unit cell.
            if (size(L) .ne. U%dim) then
                error stop "error stop in procedure new_lattice from module lattice_mod: lattice L extent dimension does not match UnitCell dimension."
            endif

            ! Check that the provided periodic extent matches the dimension of L (which matches the dimension of the unit cell).
            if (size(periodic) .ne. size(L)) then
                error stop "error stop in procedure new_lattice from module lattice_mod: lattice periodic dimension does not match lattice L extent dimension."
            endif

            Lat%U        = U
            Lat%L        = L
            Lat%periodic = periodic
            Lat%dim      = size(L)

            Lat%ncells = 1
            do i = 1, Lat%dim
                Lat%ncells = Lat%ncells * L(i)
            enddo

            Lat%nsites = Lat%ncells * U%norbitals
            Lat%nbonds = 0
            Lat%nhoppings = 0
            Lat%complexhoppings = .false.
        endfunction new_lattice

        function cind_sorb_to_sind(self, sorb, cind) result(sind)
            class(Lattice), intent(in) :: self
            integer      , intent(in) :: sorb
            integer      , intent(in) :: cind

            integer :: sind

            sind = sorb + self%U%norbitals * (cind - 1)
        endfunction cind_sorb_to_sind

        function ccds_sorb_to_sind(self, ccds, sorb) result(sind)
            class(Lattice), intent(in) :: self
            integer       , intent(in) :: ccds(self%dim)
            integer       , intent(in) :: sorb
            
            integer :: sind, cind

            cind = self%ccds_to_cind(ccds)
            sind = self%cind_sorb_to_sind(sorb, cind)
        endfunction ccds_sorb_to_sind

        function sind_to_cind(self, sind) result(cind)
            class(Lattice), intent(in) :: self
            integer      , intent(in) :: sind

            integer :: cind

            cind = (sind - 1) / self%U%norbitals+1
        endfunction sind_to_cind

        function ccds_to_cind(self, ccds) result(cind)
            class(Lattice), intent(in) :: self
            integer       , intent(in) :: ccds(:)

            integer :: cind
            integer :: i, strider
            associate(dim => self%dim, L => self%L)
                cind = 1 ; strider = 1
                do i = 1, dim
                    cind = cind + ccds(i) * strider
                    strider = strider * L(i)
                enddo
            endassociate
        endfunction ccds_to_cind

        function new_bond(orbto, orbfr, cdsp) result(B)
            integer, intent(in) :: orbto
            integer, intent(in) :: orbfr
            integer, intent(in) :: cdsp(:)

            type(Bond) :: B

            B%orbto = orbto
            B%orbfr = orbfr
            B%cdsp  = cdsp
        endfunction new_bond

        function new_unitcell(dim) result(U)
            integer, intent(in) :: dim

            type(UnitCell) :: U

            U%dim = dim
            U%norbitals = 0
            U%V = 0.0_dp
            allocate(U%lvecs(dim, dim))
            U%lvecs = 0.0_dp
        endfunction new_unitcell

        subroutine set_lvec(self, i, a)
            class(UnitCell), intent(inout) :: self
            integer        , intent(in)    :: i
            real(dp)       , intent(in)    :: a(:)

            associate(dim => self%dim, lvecs => self%lvecs)
                ! Make sure the dimension of a matches the dimension of the UnitCell.
                if (size(a) .ne. dim) stop "error stop in procedure set_lvec from module lattice_mod: mismatch in UnitCell and input lattice vector dimension."
                ! Make sure i is between 1 and dim.
                if ((i .le. 0) .or. (i .gt. dim)) error stop "error stop in procedure set_lvec from module lattice_mod: attempting to set lattice vector index out of range of UnitCell dimension."
                lvecs(:, i) = a
            endassociate
        endsubroutine set_lvec

        function lvec(self, i) result(a)
            class(UnitCell), intent(in) :: self
            integer        , intent(in) :: i

            real(dp) :: a(self%dim)

            associate(dim => self%dim, lvecs => self%lvecs)
                ! Make sure i is between 1 and dim.
                if ((i .le. 0) .or. (i .gt. dim)) error stop "error stop in procedure lvec from module lattice_mod: attempting to access lattice vector index out of range of UnitCell dimension."
                a = lvecs(:, i)
            endassociate
        endfunction lvec

        !> Get the `i`th reciprocal lattice vector of the unit cell.
        !!
        !! @param[in] self Unit cell to get reciprocal lattice vector of.
        !! @param[in] i index of reciprocal lattice vector to get.
        function rvec(self, i) result(b)
            class(UnitCell), intent(in) :: self
            integer        , intent(in) :: i

            real(dp) :: b(self%dim)

            associate(dim => self%dim, rvecs => self%rvecs)
                ! Make sure i is between 1 and dim.
                if ((i .le. 0) .or. (i .gt. dim)) error stop "error stop in procedure rvec from module lattice_mod: attempting to access lattice vector index out of range of UnitCell dimension."
                b = rvecs(:, i)
            endassociate
        endfunction rvec

        !> Adds an orbital to the unit cell at the specified position.
        !!
        !! The position vector @c r should be written in terms of the lattice vectors, with each entry `0 <= r(i) <= 1`.
        !! Within a unit cell, orbitals are added in increasing index. That is, the first orbital added will have
        !! index 1, the second index 2, and so on.
        !!
        !! @param[in,out] self Unit cell to add an orbital to.
        !! @param[in]     r position of orbital in unit cell.
        subroutine add_orbital(self, r)
            class(UnitCell), intent(inout) :: self
            real(dp)       , intent(in)    :: r(:)

            associate (dim => self%dim, norbitals => self%norbitals)
                ! Make sure the dimension of position matches the dimension of the UnitCell.
                if (size(r) .ne. dim) stop "error stop in procedure add_orbital from module lattice_mod: mismatch in UnitCell and input position vector dimension."
                call append_column(self%orbital_positions, r)
                norbitals = norbitals + 1
            endassociate
        endsubroutine add_orbital

        !> Rotates a 2-dimensional vector \f$ 90^\circ \f$ clockwise.
        !!
        !! @param[in] x vector to rotate.
        function twodim_90degree_rotation(x) result(y)
            real(dp), intent(in) :: x(2)

            real(dp) :: y(2)

            y(1) =  x(2)
            y(2) = -x(1)
        endfunction twodim_90degree_rotation

        !> Returns the cross product \f$ c = a \times b \f$.
        !!
        !! @param[in] a first vector in cross product.
        !! @param[in] b second vector in cross product.
        function cross_product(a, b) result(c)
            real(dp), intent(in) :: a(3)
            real(dp), intent(in) :: b(3)

            real(dp) :: c(3)

            c(1) = a(2)*b(3) - a(3)*b(2)
            c(2) = a(3)*b(1) - a(1)*b(3)
            c(3) = a(1)*b(2) - a(2)*b(1)
        endfunction cross_product

        subroutine make_rvecs(self)
            class(UnitCell), intent(inout) :: self

            associate (dim => self%dim, lvecs => self%lvecs, V => self%V)
                if (allocated(self%rvecs)) deallocate(self%rvecs)
                allocate(self%rvecs, mold=lvecs)

                selectcase (dim)

                    case (1)
                        associate(a1 => lvecs(:, 1), b1 => self%rvecs(:, 1))
                            V = a1(1)
                            b1 = 2.0_dp * pi * a1 / dot_product(a1, a1)
                        endassociate

                    case (2)
                        associate(a1 => lvecs(:, 1), b1 => self%rvecs(:, 1),   &
                                  a2 => lvecs(:, 2), b2 => self%rvecs(:, 2))
                            V = dot_product(a1, twodim_90degree_rotation(a2))
                            b1 = 2.0_dp * pi * twodim_90degree_rotation(a2) /   &
                                 dot_product(a1, twodim_90degree_rotation(a2))
                            b2 = 2.0_dp * pi * twodim_90degree_rotation(a1) /   &
                                 dot_product(a2, twodim_90degree_rotation(a1))
                        endassociate

                    case (3)
                        associate(a1 => lvecs(:, 1), b1 => self%rvecs(:, 1),   &
                                  a2 => lvecs(:, 2), b2 => self%rvecs(:, 2),   &
                                  a3 => lvecs(:, 3), b3 => self%rvecs(:, 3))
                            V = dot_product(a1, cross_product(a2, a3))
                            b1 = (2.0_dp * pi / V) * cross_product(a2, a3)
                            b2 = (2.0_dp * pi / V) * cross_product(a3, a1)
                            b3 = (2.0_dp * pi / V) * cross_product(a1, a2)
                        endassociate
                            
                endselect

            endassociate
        endsubroutine make_rvecs

        !> Appends a column to the @c real(dp) matrix @p A
        !!
        !! @param[in,out] A Matrix to append a column to.
        !! @param[in]     x Vector to append to @p A.
        subroutine append_column_dp(A, x)
            real(dp), allocatable, intent(inout) :: A(:, :)
            real(dp),              intent(in)    :: x(:)

            real(dp), allocatable :: temp(:, :)

            ! Make x the first column of A if A is not already allocated
            if (.not. allocated(A)) then
                allocate(A(size(x), 1))
                A(:, 1) = x
                return
            endif

            ! Make sure A and x agree in dimension.
            if (size(A, 1) .ne. size(x)) error stop "error stop in procedure append_column_dp from module lattice_mod: attempting to append a column to a matrix with mismatching number of rows."

            allocate(temp(size(x), size(A, 2) + 1))
            temp(:, 1:size(A, 2)  ) = A
            temp(:, size(A, 2) + 1) = x
            call move_alloc(temp, A)
        endsubroutine append_column_dp

        ! Fortran version of numpy's isclose
        ! https://numpy.org/doc/stable/reference/generated/numpy.isclose.html
        pure elemental function isclose_dp(a, b, rtol, atol, equal_nan) result(val)
            real(dp), intent(in)           :: a
            real(dp), intent(in)           :: b
            real(dp), intent(in), optional :: rtol
            real(dp), intent(in), optional :: atol
            logical,  intent(in), optional :: equal_nan

            logical  :: val
            real(dp) :: artol
            real(dp) :: aatol
            logical  :: aequal_nan
            logical  :: a_nan, b_nan

            ! Default values
            artol      = 1.0e-5_dp ! Relative tolerance
            aatol      = 1.0e-8_dp ! Absolute tolerance
            aequal_nan = .false.   ! Whether or not nan is treated as equal

            if (present(rtol))      artol      = rtol
            if (present(atol))      aatol      = atol
            if (present(equal_nan)) aequal_nan = equal_nan

            ! Check if input is nan
            a_nan = ieee_is_nan(a)
            b_nan = ieee_is_nan(b)
            if (a_nan .or. b_nan) then
                val = aequal_nan .and. a_nan .and. b_nan
                return
            endif

            ! Check actual equality
            if (a .eq. b) then
                val = .true.
                return
            endif

            ! Check if input is infinite
            if ((.not. ieee_is_finite(a)) .or. (.not. ieee_is_finite(b))) then
                val = .false.
                return
            end if

            ! Now guaranteed a and b are not nan or inf
            ! Check if they are within an absolute tolerance + relative tolerance.
            if (abs(a - b) .le. aatol + artol * abs(b)) then
                val = .true.
            else
                val = .false.
            endif
        end function isclose_dp

        pure elemental function isclose_zp(a, b, rtol, atol, equal_nan) result(val)
            complex(dp), intent(in)           :: a
            complex(dp), intent(in)           :: b
            real(dp)   , intent(in), optional :: rtol
            real(dp)   , intent(in), optional :: atol
            logical    , intent(in), optional :: equal_nan

            logical  :: val
            real(dp) :: artol
            real(dp) :: aatol
            logical  :: aequal_nan
            logical  :: a_nan, b_nan
            logical  :: a_finite, b_finite

            ! Default values
            artol      = 1.0e-5_dp ! Relative tolerance
            aatol      = 1.0e-8_dp ! Absolute tolerance
            aequal_nan = .false.   ! Whether or not nan is treated as equal

            if (present(rtol))      artol      = rtol
            if (present(atol))      aatol      = atol
            if (present(equal_nan)) aequal_nan = equal_nan

            ! Check if input is nan
            a_nan = ieee_is_nan(real(a, dp)) .or. ieee_is_nan(aimag(a))
            b_nan = ieee_is_nan(real(b, dp)) .or. ieee_is_nan(aimag(b))

            if (a_nan .or. b_nan) then
                val = aequal_nan .and. a_nan .and. b_nan
                return
            endif

            ! Check actual equality
            if (a .eq. b) then
                val = .true.
                return
            endif

            ! Check if input is finite
            a_finite = ieee_is_finite(real(a, dp)) .and. ieee_is_finite(aimag(a))
            b_finite = ieee_is_finite(real(b, dp)) .and. ieee_is_finite(aimag(b))

            if ((.not. a_finite) .or. (.not. b_finite)) then
                val = .false.
                return
            endif

            ! Now guaranteed a and b are not nan or inf
            ! Check if they are within an absolute tolerance + relative tolerance.
            if (abs(a - b) .le. aatol + artol * abs(b)) then
                val = .true.
            else
                val = .false.
            endif
        endfunction isclose_zp

        pure elemental function iscomplex_zp(x, rtol, atol) result(val)
            complex(dp), intent(in)           :: x
            real(dp)   , intent(in), optional :: rtol
            real(dp)   , intent(in), optional :: atol

            logical :: val

            val = .not. isclose(aimag(x), 0.0_dp, rtol=rtol, atol=atol)
        endfunction iscomplex_zp

        function new_hopping(tto, tfr) result(t)
            complex(dp), intent(in) :: tto
            complex(dp), intent(in) :: tfr

            type(Hopping) :: t

            t%tto = tto
            t%tfr = tfr
            t%iscomplex = iscomplex(tto) .or. iscomplex(tfr)
        endfunction new_hopping

        function ccds_ccds_to_cdsp(self, ccdsfr, ccdsto) result(cdsp)
            class(Lattice), intent(in) :: self
            integer       , intent(in) :: ccdsfr(self%dim)
            integer       , intent(in) :: ccdsto(self%dim)

            integer :: cdsp(self%dim)
            integer :: i

            associate (dim => self%dim, L => self%L, periodic => self%periodic)
                cdsp = ccdsto - ccdsfr
                do i = 1, dim
                    if (periodic(i)) then
                        cdsp(i) = modulo(cdsp(i), L(i))
                    endif
                enddo
            endassociate
        endfunction ccds_ccds_to_cdsp

        function cind_cind_to_cdsp(self, cindfr, cindto) result(cdsp)
            class(Lattice), intent(in) :: self
            integer       , intent(in) :: cindfr
            integer       , intent(in) :: cindto

            integer :: cdsp(self%dim)
            integer :: ccdsfr(self%dim)
            integer :: ccdsto(self%dim)

            ccdsfr = self%cind_to_ccds(cindfr)
            ccdsto = self%cind_to_ccds(cindto)
            cdsp   = self%ccds_ccds_to_cdsp(ccdsfr, ccdsto)
        endfunction cind_cind_to_cdsp

        function sind_sind_to_cdsp(self, sindfr, sindto) result(cdsp)
            class(Lattice), intent(in) :: self
            integer       , intent(in) :: sindfr
            integer       , intent(in) :: sindto

            integer :: cdsp(self%dim)
            integer :: cindfr
            integer :: cindto

            cindfr = self%sind_to_cind(sindfr)
            cindto = self%sind_to_cind(sindto)
            cdsp = self%cind_cind_to_cdsp(cindfr, cindto)
        endfunction sind_sind_to_cdsp

        function oind_cdsp_oind_to_sep(self, oindfr, cdsp, oindto) result(sep)
            class(Lattice), intent(in) :: self
            integer       , intent(in) :: oindfr
            integer       , intent(in) :: cdsp(self%dim)
            integer       , intent(in) :: oindto

            real(dp) :: sep(self%dim)
            real(dp) :: sdsp(self%dim)

            sdsp = real(cdsp, dp) + self%U%orbital_positions(:, oindto) - self%U%orbital_positions(:, oindfr)
            sep = self%cds_to_pos(sdsp)
        endfunction oind_cdsp_oind_to_sep

        function sind_sind_to_sep(self, sindfr, sindto) result(sep)
            class(Lattice), intent(in) :: self
            integer       , intent(in) :: sindfr
            integer       , intent(in) :: sindto

            real(dp) :: sep(self%dim)

            integer :: cindfr, cindto
            integer :: oindfr, oindto
            integer :: cdsp(self%dim)

            call self%sind_to_cind_oind(sindfr, cindfr, oindfr)
            call self%sind_to_cind_oind(sindto, cindto, oindto)

            cdsp = self%cind_cind_to_cdsp(cindfr, cindto)

            sep = self%oind_cdsp_oind_to_sep(oindfr, cdsp, oindto)
        endfunction sind_sind_to_sep

        subroutine complexhoppings_check(self)
            class(Lattice), intent(inout) :: self

            integer :: i

            self%complexhoppings = .false.
            do i = 1, self%nhoppings
                if (self%hoppings(i)%iscomplex) then
                    self%complexhoppings = .true.
                    return
                endif
            enddo
        endsubroutine complexhoppings_check






        function hopping_matrix_zp(self) result(K)
            class(Lattice), intent(in) :: self

            complex(dp), allocatable :: K(:, :)
            type(Bond)    :: B
            type(Hopping) :: t

            integer :: bind
            integer :: cindfr
            integer :: sindfr
            integer :: sindto

            logical :: in_lattice

            allocate(K(self%nsites, self%nsites))
            K = 0.0_dp

            ! Loop over all bonds
            do bind = 1, self%nbonds

                B = self%bonds(bind)
                t = self%hoppings(bind)

                ! Loop over all cells 
                do cindfr = 1, self%ncells

                    ! Get the site index the current bond is from in the current cell.
                    sindfr = self%cind_sorb_to_sind(B%orbfr, cindfr)

                    ! Get the site index the bond goes to.
                    call self%cind_cdsp_sorb_to_sind(cindfr, B%cdsp, B%orbto, sindto, in_lattice)
                    
                    ! Make sure the last result was inside of the lattice.
                    if (.not. in_lattice) cycle

                    ! fr --> to.
                    K(sindto, sindfr) = K(sindto, sindfr) + t%tto

                    ! fr <-- to.
                    K(sindfr, sindto) = K(sindfr, sindto) + t%tfr

                enddo
            enddo
        endfunction hopping_matrix_zp







endmodule lattice_mod