
c     The number of shells
      integer NSHELL
c     The number of primitives (total)
      integer NPRIMT
c     The number of primitives squared (total)
      integer NPRIM2
c     The number of unique centers
      integer NUCENT

c     The pointer to an array with the exponents (nprimt)
      integer PEXPON
c     The pointer to an array with the contraction coeffs (nprimt)
      integer PCONCO
c     The pointer to an array with the types (nprimt)
      integer PTYPES
c     The pointer to where the shell for this center begins (nucent)
      integer PPSHC
c     The pointer to an array which gives the number of shells on each center
      integer PNSHC
c     The pointer to an array with the number of prims in each shell
      integer PNPRIM
c     Pointer to an array with pointers to the list of primitives for each shell
      integer PPPRIM
c     Pointer to an array with pointers to the list of primitives for each shell
      integer PCHARG
c     Pointer to an array with the geometry. Dimensioned as geom(3,nucent)
      integer PGEOM

      common /intbs1/ nshell,nprimt,nprim2,nucent
     &   , pexpon, pconco, ptypes, ppshc, pnshc, pnprim, ppprim
     &   , pcharg, pgeom

