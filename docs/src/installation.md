# Installation

TPLib.jl is a binding to the library [TPLib](http://www.cmap.polytechnique.fr/~allamigeon/software/) written in OCaml. You must therefore install the library on your system. You can download it here and build it from source:

```
./configure
make
make install
```

You can also install it with the OCaml package manager [OPAM](https://opam.ocaml.org):

```
opam install tplib
```

Computing with `Rationals` or `BigInts` requires the [zarith](https://opam.ocaml.org/packages/zarith/) package for OCaml. It can be installed using OPAM:

```
opam install zarith
```
