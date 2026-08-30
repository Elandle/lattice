!
! gfortran -c -ffree-line-length-none lattice_mod.f90
!
module lattice_mod
    use iso_fortran_env, only: dp => real64
    implicit none

    real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)

    type :: UnitCell
        integer  :: dim       ! Dimension of unit cell
        integer  :: norbitals ! Number of orbitals in the unit cell
        real(dp) :: V         ! Volume 

        real(dp), allocatable :: lattice_vectors(:, :)
        real(dp), allocatable :: orbital_positions(:, :)
        real(dp), allocatable :: reciprocal_vectors(:, :)

        contains
            procedure set_lattice_vector
            procedure lattice_vector
            procedure add_orbital
            procedure make_reciprocal_vectors
            procedure reciprocal_vector
    endtype UnitCell

    type :: Bond
        integer :: to   ! Orbital index the bond is to.
        integer :: from ! Orbital index the bond is from.   

        integer, allocatable :: displacement(:) ! Unit cell displacement.
                                                ! Measured in terms of lattice vectors.
    endtype Bond

    type :: Hopping
        type(Bond)  :: B
        complex(dp) :: t_to   ! Hopping from  ---  t_to   ---> to
        complex(dp) :: t_from ! Hopping from <--- t_from  ---  to

        ! Internally store hoppings as complex numbers, then check later on if they are real or complex
        logical :: iscomplex
    endtype Hopping


    type :: Lattice
        type(UnitCell)              :: U           ! Unit cell of the lattice.
        integer       , allocatable :: L(:)        ! Integer extent of lattice in units of primitive vectors.
                                                   ! L(i) = number of unit cells along axis defined by lattice vector i.
        logical       , allocatable :: periodic(:) ! Which primitive vector directions are periodic.
        integer                     :: ncells      ! Number of unit cells in the lattice.
                                                   ! ncells = product of L entries (eg, 2 dimensional ncells = L(1) * L(2))
        integer                     :: nsites      ! Number of sites. nsites = ncells * norbitals (norbitals = number of orbitals per site).
        integer                     :: dim
        type(Bond)    , allocatable :: bonds(:)
        integer                     :: nbonds

    contains
        procedure :: siteindx_from_cellindx
        procedure :: siteindx_from_cellcoords
        procedure :: cellindx_from_siteindx
        procedure :: cellindx_from_coords
        procedure :: cell_coordinatess_from_indx

        procedure :: cell_indx_displacement
        procedure :: cell_displacement => cell_coordinate_displacement

        procedure :: add_bond
        procedure :: site_indx_displacement
    endtype Lattice

    interface append_column
        module procedure :: append_column_dp
    endinterface append_column

    contains

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

        subroutine site_indx_displacement(self, indxfrom, dr, orbitalto, indxto, in_lattice)
            class(Lattice), intent(in)  :: self
            integer       , intent(in)  :: indxfrom
            integer       , intent(in)  :: dr(self%dim)
            integer       , intent(in)  :: orbitalto
            integer       , intent(out) :: indxto
            logical       , intent(out) :: in_lattice

            integer :: cellfrom, cellto

            cellfrom = self%cellindx_from_siteindx(indxfrom)
            call self%cell_indx_displacement(cellfrom, dr, cellto, in_lattice)
            if (in_lattice) indxto = self%siteindx_from_cellindx(orbitalto, cellto)
        endsubroutine site_indx_displacement

        subroutine cell_indx_displacement(self, indxfrom, dr, indxto, in_lattice)
            class(Lattice), intent(in)  :: self
            integer       , intent(in)  :: indxfrom
            integer       , intent(in)  :: dr(self%dim)
            integer       , intent(out) :: indxto
            logical       , intent(out) :: in_lattice

                integer :: rfrom(self%dim)
                integer :: rto(self%dim)

                rfrom = self%cell_coordinatess_from_indx(indxfrom)

                call self%cell_displacement(rfrom, dr, rto, in_lattice)

                if (in_lattice) indxto = self%cellindx_from_coords(rto)
        endsubroutine cell_indx_displacement

        subroutine cell_coordinate_displacement(self, rfrom, dr, rto, in_lattice)
            ! Sets rto = rfrom + dr in lattice coordinates keeping track of periodicity.
            ! If rto is out of bounds for a non-periodic direction, in_lattice is set to .false.
            ! otherwise (rto is inside of the lattice) in_lattice is set to .true. .
            class(Lattice), intent(in)  :: self
            integer       , intent(in)  :: rfrom(self%dim)
            integer       , intent(in)  :: dr(self%dim)
            integer       , intent(out) :: rto(self%dim)
            logical       , intent(out) :: in_lattice

            integer :: i

            associate (dim => self%dim, L => self%L, periodic => self%periodic)

                rto = rfrom + dr
                in_lattice = .true.

                do i = 1, dim
                    if (periodic(i)) then
                        rto(i) = modulo(rto(i), L(i))
                    elseif ((rto(i) .lt. 0) .or. (rto(i) .ge. L(i))) then
                        in_lattice = .false.
                        return
                    endif
                enddo
            endassociate
        endsubroutine cell_coordinate_displacement

        function cell_coordinatess_from_indx(self, cellindx) result(r)
            class(Lattice), intent(in) :: self
            integer, intent(in) :: cellindx

            integer :: r(self%dim)
            integer :: i, n

            n = cellindx - 1

            do i = 1, self%dim
                r(i) = modulo(n, self%L(i))
                n = n / self%L(i)
            enddo
        endfunction cell_coordinatess_from_indx




        function new_lattice(U, L, periodic) result(Lat)
            type(UnitCell), intent(in) :: U
            integer       , intent(in) :: L(:)
            logical       , intent(in) :: periodic(:)

            type(Lattice) :: Lat

            integer :: i


            Lat%U = U
            Lat%L = L
            Lat%periodic = periodic
            Lat%dim = size(L)

            Lat%ncells = 1
            do i = 1, Lat%dim
                Lat%ncells = Lat%ncells * L(i)
            enddo
            Lat%nsites = Lat%ncells * U%norbitals
        endfunction new_lattice

        function siteindx_from_cellindx(self, orbital, cellindx) result(siteindx)
            ! Index of a site given its orbital number and the index of the unit cell it is in.
            class(Lattice), intent(in) :: self
            integer      , intent(in) :: orbital
            integer      , intent(in) :: cellindx

            integer :: siteindx

            siteindx = orbital + self%U%norbitals * (cellindx - 1)
        endfunction siteindx_from_cellindx

        function siteindx_from_cellcoords(self, orbital, r) result(siteindx)
            ! Index of a site given its orbital number and the coordinates of the unit cell it is in.
            class(Lattice), intent(in) :: self
            integer      , intent(in) :: orbital
            integer      , intent(in) :: r(:)

            integer :: siteindx
            integer :: cellindx

            cellindx = self%cellindx_from_coords(r)
            siteindx = self%siteindx_from_cellindx(orbital, cellindx)
        endfunction siteindx_from_cellcoords

        function cellindx_from_siteindx(self, siteindx) result(cellindx)
            ! Index of the unit cell a site is in.
            class(Lattice), intent(in) :: self
            integer      , intent(in) :: siteindx

            integer :: cellindx

            cellindx = (siteindx - 1) / self%U%norbitals+1
        endfunction cellindx_from_siteindx

        function cellindx_from_coords(self, r) result(cellindx)
            ! Index of a unit cell at a specified coordinate.
            class(Lattice), intent(in) :: self
            integer      , intent(in) :: r(:)

            integer :: cellindx
            integer :: i, strider
            associate(dim => self%dim, L => self%L)
                cellindx = 1 ; strider = 1
                do i = 1, dim
                    cellindx = cellindx + r(i) * strider
                    strider = strider * L(i)
                enddo
            endassociate
        endfunction cellindx_from_coords


        function new_bond(to, from, displacement) result(B)
            integer, intent(in) :: to
            integer, intent(in) :: from
            integer, intent(in) :: displacement(:)

            type(Bond) :: B

            B%to = to
            B%from = from
            B%displacement = displacement
        endfunction new_bond

        function new_unitcell(dim) result(U)
            integer, intent(in) :: dim

            type(UnitCell) :: U

            U%dim = dim
            U%norbitals = 0
            U%V = 0.0_dp
            allocate(U%lattice_vectors(dim, dim))
            U%lattice_vectors = 0.0_dp
        endfunction new_unitcell

        subroutine set_lattice_vector(self, i, a)
            class(UnitCell), intent(inout) :: self
            integer        , intent(in)    :: i
            real(dp)       , intent(in)    :: a(:)

            associate(dim => self%dim, lattice_vectors => self%lattice_vectors)
                ! Make sure the dimension of a matches the dimension of the UnitCell.
                if (size(a) .ne. dim) stop "error stop in procedure set_lattice_vector from module lattice_mod: mismatch in UnitCell and input lattice vector dimension."
                ! Make sure i is between 1 and dim.
                if ((i .le. 0) .or. (i .gt. dim)) error stop "error stop in procedure set_lattice_vector from module lattice_mod: attempting to set lattice vector index out of range of UnitCell dimension."
                lattice_vectors(:, i) = a
            endassociate
        endsubroutine set_lattice_vector

        function lattice_vector(self, i) result(a)
            class(UnitCell), intent(in) :: self
            integer        , intent(in) :: i

            real(dp) :: a(self%dim)

            associate(dim => self%dim, lattice_vectors => self%lattice_vectors)
                ! Make sure i is between 1 and dim.
                if ((i .le. 0) .or. (i .gt. dim)) error stop "error stop in procedure lattice_vector from module lattice_mod: attempting to access lattice vector index out of range of UnitCell dimension."
                a = lattice_vectors(:, i)
            endassociate
        endfunction lattice_vector

        !> Get the `i`th reciprocal lattice vector of the unit cell.
        !!
        !! @param[in] self Unit cell to get reciprocal lattice vector of.
        !! @param[in] i index of reciprocal lattice vector to get.
        function reciprocal_vector(self, i) result(b)
            class(UnitCell), intent(in) :: self
            integer        , intent(in) :: i

            real(dp) :: b(self%dim)

            associate(dim => self%dim, reciprocal_vectors => self%reciprocal_vectors)
                ! Make sure i is between 1 and dim.
                if ((i .le. 0) .or. (i .gt. dim)) error stop "error stop in procedure reciprocal_vector from module lattice_mod: attempting to access lattice vector index out of range of UnitCell dimension."
                b = reciprocal_vectors(:, i)
            endassociate
        endfunction reciprocal_vector

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

        subroutine make_reciprocal_vectors(self)
            class(UnitCell), intent(inout) :: self

            associate (dim => self%dim, lattice_vectors => self%lattice_vectors, V => self%V)
                if (allocated(self%reciprocal_vectors)) deallocate(self%reciprocal_vectors)
                allocate(self%reciprocal_vectors, mold=lattice_vectors)

                selectcase (dim)

                    case (1)
                        associate(a1 => lattice_vectors(:, 1), b1 => self%reciprocal_vectors(:, 1))
                            V = a1(1)
                            b1 = 2.0_dp * pi * a1 / dot_product(a1, a1)
                        endassociate

                    case (2)
                        associate(a1 => lattice_vectors(:, 1), b1 => self%reciprocal_vectors(:, 1),   &
                                  a2 => lattice_vectors(:, 2), b2 => self%reciprocal_vectors(:, 2))
                            V = dot_product(a1, twodim_90degree_rotation(a2))
                            b1 = 2.0_dp * pi * twodim_90degree_rotation(a2) /   &
                                 dot_product(a1, twodim_90degree_rotation(a2))
                            b2 = 2.0_dp * pi * twodim_90degree_rotation(a1) /   &
                                 dot_product(a2, twodim_90degree_rotation(a1))
                        endassociate

                    case (3)
                        associate(a1 => lattice_vectors(:, 1), b1 => self%reciprocal_vectors(:, 1),   &
                                  a2 => lattice_vectors(:, 2), b2 => self%reciprocal_vectors(:, 2),   &
                                  a3 => lattice_vectors(:, 3), b3 => self%reciprocal_vectors(:, 3))
                            V = dot_product(a1, cross_product(a2, a3))
                            b1 = (2.0_dp * pi / V) * cross_product(a2, a3)
                            b2 = (2.0_dp * pi / V) * cross_product(a3, a1)
                            b3 = (2.0_dp * pi / V) * cross_product(a1, a2)
                        endassociate
                            
                endselect

            endassociate
        endsubroutine make_reciprocal_vectors

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
endmodule lattice_mod