# Conventions

This page describes some of the conventions used in this library.

# Lattices, Orbitals, Bonds, and Hoppings
At the end of the day, a lattice is a collection of sites indexed by \f$ i = 1, 2, \dots, m \f$ with a notion of hoppings \f$ t_{ij} \f$ from site \f$ j \f$ to site \f$ i \f$ (eg, in the Ising model \f$ J_{ij} = t_{ij} \f$ is interpreted as an interaction strength between sites \f$ j \f$ and \f$ i \f$).

## Lattices
Lattices are usually regular in their structure, so lattice vectors \f$ \mathbf{a}_i \f$ \f$ i = 1, \dots, d \f$ are typically used to describe them (\f$ d = \f$ the dimension of the lattice, and is \f$ 1 \f$, \f$ 2 \f$, or \f$ 3 \f$ for this code). The idea is that the lattice vectors define a primitive unit cell and this primitive unit cell is duplicated \f$ L_1\dots L_d \f$ times (\f$ L_i \f$ is called the length along the \f$ i \f$th or \f$ \mathbf{a}_i \f$ axis). In this code a primitive unit cell is defined to be the space enclosed by the vectors \f$ \mathbf{a}_1, \dots, \mathbf{a}_d \f$ in \f$ d \f$-dimensional space. The origin of this primitive unit cell is taken to be the vector \f$ \mathbf{0} \f$ (sometimes interpreted as the "bottom-left corner"). Duplication of this primitive unit cell creates the unit cells of the lattice (unit cell = copy of primitive unit cell at a specific location).

The duplication process is this.
The \f$ L_1, \dots, L_d \f$ define the size of the lattice (eg, you might hear of an \f$ m \times n \f$ square lattice. Here \f$ d = 2\f$, \f$ L_1 = m \f$, \f$ L_2 = n \f$, and the unit cell is a square, typically \f$ \mathbf{a}_1 = (1, 0) \f$ and \f$ \mathbf{a}_2 = (0, 1) \f$).
The lattice consists of \f$ L_1\dots L_d \f$ unit cells each with their origin located at \f$ \mathbf{R} = n_1\mathbf{a}_1 + \dots + n_d\mathbf{a}_d \f$ for each combination of \f$ 0\leq n_i < L_i \f$.
In the code each unit cell is associated with a unique index \f$1\leq j \leq L_1\dots L_d\f$. The precise index mapping is:
\f[
j = 1 + \sum_{i=1}^{d} n_i \prod_{k=1}^{i-1} L_k = 1 + n_1 + n_2L_1 + n_3L_1L_2
\f]
where the sum on the right is cutoff based on the dimension \f$ d \f$.
This allows unit cells to be described either in terms of their index \f$ j \f$ or coordinates \f$ \mathbf{n}_j = (n_1, n_2, n_3) \f$:
\f[
    j = 1 + n_1 + n_2L_1 + n_3L_1L_2 \leftrightarrow (n_1, n_2, n_3)
\f]
(as always, this might be cutoff by the dimension \f$ d \f$).
In words: the \f$ (0, 0, 0) \f$ unit cell is assigned index \f$ 1 \f$ and indices increase with the first coordinate first, then second coordinate, and lastly the third coordinate.
Unit cell coordinates correspond to the physical positions of the unit cells. If the origin of the \f$ 1 \f$st unit cell is taken as the origin \f$ \mathbf{0} \f$ of space, then the position of the \f$ j \f$th unit cell's origin (taken to be called the position of the unit cell) is:
\f[
    \mathbf{R}_j = n_1\mathbf{a}_1 + n_2\mathbf{a}_2 + n_3\mathbf{a}_3
\f]
## Sites
Sites are incorporated by assigning what are called orbitals to unit cells. An orbital is simply a site inside (or on the boundary of) a unit cell. Orbitals consist of two parts: their index \f$ k \f$ in the unit cell and their position \f$ \mathbf{r}_k = (r_{1k}, r_{2k}, r_{3k}) \f$ within the unit cell in the lattice vector basis \f$ \mathbf{a}_i \f$. For example, a position \f$ (0, 0, 0) \f$ corresponds to an orbital at the origin or "bottom left" of a unit cell and \f$ (1, 1, 1) \f$ corresponds to the "top right" of a unit cell. Since orbital position vectors must lie in or on the boundary of unit cells and use coordinates in the lattice vector basis, their components must lie \f$ 0 \leq r_{ik} \leq 1 \f$. The physical position \f$ \mathbf{R}_{jk} \f$ of the \f$ k \f$th orbital in the \f$ j \f$th unit cell in the lattice is:
\f[
    \mathbf{R}_{jk} = \mathbf{R}_j + r_{1k}\mathbf{a}_1 + r_{2k}\mathbf{a}_2 + r_{3k}\mathbf{a}_3
\f]
where \f$ \mathbf{R}_j \f$ is the physical position of the \f$ j \f$th unit cell.
## Bonds
Unit cells and orbitals are sufficient for describing the locations of all sites in a lattice. Bonds indicate hoppings between sites. Like how unit cells were used as an easy tool to create lattices since lattices are typically regular, bonds are an easy tool to create hoppings for regular lattices. Three pieces of information describe a bond: the orbital index \f$ j \f$ a bond comes from, the orbital index \f$ i \f$ a bond goes to, and the displacement vector \f$ \mathbf{d} = (d_1, d_2, d_3) \f$ from the unit cell orbital \f$ j \f$ is in to the one orbital \f$ i \f$ is in in terms of the lattice vectors. An intention of this library is to allow asymmetric hoppings \f$ t_{ij} \neq t_{ji} \f$ so this order does matter, bonds are directional: from \f$ \to \f$ to. However, for hopping implementation (to be covered soon), only \f$ 1 \f$ bond per pair of sites should be defined. When a lattice is informed of a bond, it will form a bond between all \f$ j \f$, \f$ i \f$-indexed orbitals between unit cells a displacement \f$ \Delta R = d_1\mathbf{a}_1 + d_2\mathbf{a}_2 + d_3\mathbf{a}_3 \f$ apart. To be more specific about this displacement, take the \f$ l \f$th unit cell. It has a position \f$ \mathbf{R}_l \f$. Form the sum \f$ \mathbf{R}_l + \Delta R \f$. As long as this does not happen to hit \f$ \mathbf{R}_l \f$ or is outside of the lattice in the case of open boundary conditions (boundary conditions will be covered soon), it will hit the \f$ p \f$th unit cell with position \f$ \mathbf{R}_p = \mathbf{R}_l + \mathbf{d} \f$. It is then said that there is a bond from the \f$ j \f$th orbital in unit cell \f$ l \f$ to the \f$ i \f$th orbital in unit cell \f$ p \f$.

## Boundary conditions
Boundary conditions indicate whether or not a displacement vector \f$ \mathbf{d} \f$ can "wrap around" to the other side of the lattice. A boundary is either open or periodic along each \f$ \mathbf{a}_i \f$ axis, with periodic being more typical. Boundary conditions are best understood in terms of the lattice vector basis. Periodic boundary conditions identify:
\f\[
    \mathbf{n} = (n_1, n_2, n_3) \sim (n_1 + m_1L_1, n_2 + m_2L_2, n_3 + m_3L_3)
\f\]
where the \f$ m_i \f$ are integers.
Here periodic boundary conditions are written assuming periodicity in all lattice vector axes. Open boundary conditions do not have this identification. That is, for open boundary conditions in lattice vector axes:
\f\[
    \mathbf{n} = (n_1, n_2, n_3) \not\sim (n_1 + m_1L_1, n_2 + m_2L_2, n_3 + m_3L_3)
\f]
In this library open or periodic boundary conditions are imposed on each lattice vector axis independently. For example, the \f$ \mathbf{a}_1 \f$ and \f$ \mathbf{a}_2 \f$ can be chosen to be periodic, but \f$ \mathbf{a}_3 \f$ not. This means that the \f$ \mathbf{a}_1 \f$ and \f$ \mathbf{a}_2 \f$ axes have the periodicity identification, but not the \f$ \mathbf{a}_3 \f$ axis. For open boundary conditions, displacement might not yield a position on the lattice.

In particular for bonding, this means that it is possible for orbitals in some unit cells to not have a bond when one might be indicated for a bond. This occurs when a bond's displacement vector added to the unit cell vector of the orbital does not reach a different unit cell.

Boundary conditions are important for reciprocal lattice vectors (ie, \f$ \mathbf{k} \f$-space) and this will be discussed later.

## Hoppings
Typically hoppings are thought of as being attached to bonds, but in this library bonds are thought of as being attached to hoppings. It should be emphasized that this is a design choice of the library (and not something that is standard). A hopping serves to indicate the hopping values between sites in a particular bond. A hopping consists of three parts: a bond from an orbital \f$ j \f$ to an orbital \f$ i \f$, the hopping \f$ t_{ij} \f$ from \f$ j \f$ to \f$ i \f$, and the hopping \f$ t_{ji} \f$ from \f$ i \f$ to \f$ j \f$.

# Reciprocal Lattice
A reciprocal lattice can be defined for any lattice, but it is only useful for the periodic directions of a lattice.
For a lattice with lattice vectors \f$ \mathbf{i} \f$ for \f$ i = 1, \dots, d \f$, there are \f$ d \f$ reciprocal lattice vectors \f$ \mathbf{b}_j \f$ (withstanding open or periodic boundary conditions) defined by:
\f\[
    \mathbf{a}_i \cdot \mathbf{b}_j = 2\pi\delta_{ij}
\f]
Formulas exist for these \f$ \mathbf{b}_j \f$ vectors.
In \f$ d = 1 \f$:
\f\[
    \mathbf{b}_1 = \frac{2\pi}{V^2}\mathbf{a}_1 \\
    V = a_1
\f]
\f$ d = 2 \f$:
\f\[
    \mathbf{b}_1 = \frac{2\pi}{V}R\mathbf{a}_2 \\
    \mathbf{b}_2 = -\frac{2\pi}{V}R\mathbf{a}_1 \\
    V = \mathbf{a}_1\cdot R\mathbf{a}_2
\f]
where \f$ R \f$ is the \f$ \pi/2 \f$ clockwise rotation matrix in two dimensions. For \f$ d = 3 \f$:
\f\[
    \mathbf{b}_1 = \frac{2\pi}{V}\mathbf{a}_2\times\mathbf{a}_3 \\
    \mathbf{b}_2 = \frac{2\pi}{V}\mathbf{a}_3\times\mathbf{a}_1 \\
    \mathbf{b}_3 = \frac{2\pi}{V}\mathbf{a}_1\times\mathbf{a}_2 \\
    V = \mathbf{a}_1 \cdot (\mathbf{a}_2\times \mathbf{a}_3)
\f]
In all three cases \f$ V \f$ is the signed volume of the unit cell defined by the \f$ \mathbf{a}_i \f$.

