# Non-uniqueness for the complex Ginzburg-Landau equation

This repository contains the code for the computer-assisted parts of
the proofs for the paper [Non-uniqueness for the complex
Ginzburg-Landau equation](TODO).

The results of the paper are presented in two different
[Pluto.jl](https://plutojl.org/) notebooks, found in the
[`proof`](proof) directory. These notebooks are responsible for
generating all the numbers and figures that appear in the paper. It is
possible to view the results of notebooks without running any code by
opening the corresponding html-files found in the [`proof`](proof)
directory, they can be opened in any browser such as Firefox. The two
notebooks are:

- [unstable-eigenvalue.jl](proof/unstable-eigenvalue.jl) - Contains
  the proof of Theorems 4.1, 4.2 and 4.3, corresponding to the
  existence of a backward self-similar solution, an associated forward
  self-similar solution and an unstable eigenvalue.

## Reproducing the proofs

The proofs were generated with Julia version 1.12.6. This repository
contains the same `Manifest.toml` file as was used when running the
proofs; this allows installing exactly the same versions of the Julia
packages. Part of the computations are done using the C++ library
[CAPD](https://github.com/CAPDGroup/CAPD). There is no Julia package
wrapping the CAPD library, and programs therefore have to be compiled
outside of Julia, which complicates the setup process. Information
about how to install CAPD and compile the required code is found in
[`CAPD/README.md`](CAPD/README.md).

Once CAPD has been installed and the required programs compiled, the
Julia part of the code can be set up by starting Julia from this
directory and running

``` julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.test()
```

This will likely take some time the first time you run it since it has
to download and compile all the packages. If the tests run
successfully then you should be good to go!

To run any of the notebooks, you first need to start Pluto. Start
Julia from this directory and run:

``` julia
using Pkg
Pkg.activate(".")
using Pluto
Pluto.run()
```

which should open a Pluto tab in your browser. Now you can open the
notebooks inside the `proof` directory through this and it should
allow you to run the proof.

## Notes about implementation

The code in this repository is spread out over four directories:

1. `src/` - This directory contains the vast majority of the code.
   More information about it is given below.
2. `CAPD/src/` - This directory contains the CAPD code, more
   information about it is found in
   [`CAPD/README.md`](CAPD/README.md).
3. `proof/` - This directory contains the notebooks producing the
   final proofs.

The implementation is based on the code from
[CGL.jl](https://github.com/Joel-Dahne/CGL.jl), the repository with
the code associated with the proof of [Self-similar singular solutions
to the nonlinear Schrödinger and the complex Ginzburg-Landau
equation](https://arxiv.org/abs/2601.16285). There has, however, been
a lot of code both added and removed.

Some good information to have:
- Some of the parameters are fixed throughout the computations; these
  are $d$, $\omega$, $\sigma$ and $\delta$. To simplify the code, we
  store all of these in a struct, `CGLParams`, which is passed around
  under the name `Λ`. For this paper we always have $d = 3$, $\omega =
  1$, $\sigma = 1$ and $\delta = 0$, but we keep them as variables to
  make the code more flexible.
- Many of the implemented functions contain documentation that
  explains what they do. These occasionally refer to equations or
  lemmas in the paper. For these references we use the notation "REF",
  e.g., "Lemma REF(lemma:U)". We use labels, e.g., `lemma:U`, rather
  than numbers to make it easier to keep the code and the paper in
  sync. The labels are not directly seen in the PDF version of the
  paper, but can be accessed by downloading the LaTeX code for the
  paper from ArXiv. **TODO:** Discuss labels for old vs new paper.
- Many of the functions have associated tests in [`test`](test) that
  serve to increase the confidence in the implementation. All of these
  tests can be run with `Pkg.test()` as discussed above.

Section 3 of the paper how proving existence of forward and backward
self-similar solutions as well as an unstable eigenvalue can be
reduced to proving the existence of three different functions,

$$G(\mu, \gamma, \kappa),\quad \hat{G}(\nu, \hat{\gamma}_{2})\quad \text{and} H(\lambda).$$

These are implemented in the three files:

- [`G.jl`](src/G.jl)
- [`G_hat.jl`](src/G_hat.jl)
- [`H.jl`](src/H.jl)

For $G$ we need to compute enclosures of $Q_0$ and $Q_\infty$, for
$\hat{G}$ enclosures of $\hat{Q}_0$ and $\hat{Q}_\infty$ and for $H$ of
$Y_0$ and $Y_\infty$. Again, see Section 3 for more details. The
implementations of these functions is what encompasses most of the
code in this repository and is split over six directories:

- [`Q_zero/`](src/Q_zero) - Related to Section 8 in the self-similar
  blowup paper. Consists of the 5 files:
  - `equation.jl` - Implements the equation for the ODE as
     well as the recursion for computing the coefficients of the
     Taylor expansion.
  - `Q_float.jl` - Contains the implementation of the non-rigorous
     evaluation of $Q_0$.
  - `Q_zero/Q_taylor.jl` - Contains the code for enclosing $Q_0$ using
     the Taylor expansion at zero.
  - `Q_zero/Q_capd.jl` - Contains the Julia code that is responsible
     for calling the CAPD code, which then computes enclosures of
     $Q_0$.
  - `Q_zero/Q.jl` - Defines wrapper methods that call either the
     `Q_float.jl` methods or `Q_capd.jl` methods depending on the type
     of the input. This is the interface that is primarily used by the
     other parts of the code.
- [`Q_infinity/`](src/Q_infinity) - Related to Section 7 in the
  self-similar blowup paper. Consists of the 9 files:
  - `parameters.jl` - Contains code for computing $a$, $b$, $c$,
    $B_W$, $B_{W,\kappa}$ and $B_{W,\epsilon}$.
  - `functions.jl` - Contains code for evaluating $P(\xi)$, $E(\xi)$
    and related functions. The functions follow the same naming scheme
    as in the paper, with the derivatives denoted by a postfix, for
    example, the function `P_dξ_dξ_dκ` is used to compute
    $P_\kappa''$.
  - `function_bounds.jl` - Responsible for computing the bounds that
    occur in Lemma 7.3 in the paper.
  - `I_bounds.jl` - The bounds from Sections 7.3.2 and 7.3.3 in the
    paper are implemented here.
  - `norm_bounds.jl` and `norm_bounds_constants.jl` - Contain the code
    for computing initial bounds of the norms of $Q_\infty$ and its
    derivatives.
  - `I_P.jl` - This code computes enclosures, not only bounds, for
    $I_P$, $I_{P,\gamma}$, $I_{P,\kappa}$, and $I_{P,\epsilon}$. It is
    based on the approach discussed in Section 7.3.4.
  - `Q.jl` - Puts all of the above code together to compute enclosures
    for $Q$ and its Jacobian. It defines the functions `Q_infinity`,
    `Q_infinity_jacobian_kappa` and `Q_infinity_jacobian_epsilon`.
- [`Q_hat_zero/`](src/Q_hat_zero) - Related to Section 6 of the paper.
  To a large part it reuses the implementation in `Q_zero/`. It
  consists of the 4 files:
  - `Q_hat_float.jl` - Contains the implementation of the non-rigorous
     evaluation of $\hat{Q}_0$.
  - `Q_hat_taylor.jl` - Contains the code for enclosing $\hat{Q}_0$ using
     the Taylor expansion at zero.
  - `Q_hat_capd.jl` - Contains the Julia code that is responsible for
     calling the CAPD code, which then computes enclosures of
     $\hat{Q}_0$.
  - `Q_hat.jl` - Defines wrapper methods that call either the
     `Q_hat_float.jl` methods or `Q_hat_capd.jl` methods depending on
     the type of the input. This is the interface that is primarily
     used by the other parts of the code.
- [`Q_hat_infinity/`](src/Q_hat_infinity) - Related to Section 5 of
  the paper. Consists of the 6 files:
  - `functions.jl` - Contains code for evaluating $\hat{P}(\xi)$,
    $\hat{E}(\xi)$ and related functions.
  - `function_bounds.jl` - Responsible for computing the bounds that
    occur in Lemma 5.1, 5.2 and 5.3 in the paper.
  - `norm_bounds.jl` - Contain the code for computing initial bounds
    of the norms of $\hat{Q}_\infty$ and $\hat{Q}_\infty'$.
  - `Q_hat_bounds.jl` - Computes bounds related to the leading
    asymptotic behavior of $\hat{Q}_\infty$ from Lemma 5.5.
  - `I_E_hat.jl` - Computes enclosures, not only bounds, for
    $I_\hat{E}$, and its derivatives with respect to the real and
    imaginary parts of `\hat{\gamma}_2`. It is based on Lemma 5.6 and
    5.7.
  - `Q_hat.jl` - Puts all of the above code together to compute
    enclosures for $\hat{Q}$ and its Jacobian.
- [`Y_zero/`](src/Y_zero) - Related to Section 8 of the paper.
  Consists of the 5 files:
  - `equation.jl` - TODO
  - `Y_float.jl` - TODO
  - `Y_taylor.jl` - TODO
  - `Y_capd.jl` - TODO
  - `Y.jl` - TODO
- [`Y_infinity/`](src/Y_infinity) - Related to Section 7 of the paper.
  Consists of the 5 files:
  - `functions.jl` - TODO
  - `function_bounds.jl` - TODO
  - `norm_bounds.jl` - TODO
  - `I_K_2.jl` - TODO
  - `Y.jl` - TODO

The implementations in `Q_capd.jl`, `Q_hat_capd.jl` and `Y_capd.jl`
make use of a rigorous integration implemented using the CAPD library.
This is implemented in two files:

- [`Q_zero.cpp`](CAPD/src/Q_zero.cpp) - Handles enclosures of both $Q_0$
  and $\hat{Q}_0$.
- [`Y_zero.cpp`](CAPD/src/Y_zero.cpp) - Handles enclosures of $Y_0$.

Apart from the above mentioned files and directories, there are a
number of auxiliary files, most in the [`src/`](src) directory,
handling various tasks:

- `assert_proof.jl`, `arb.jl`, `helper.jl`, `special-functions.jl` -
  These files implement a variety of convenience functions, such as
  formatting of output and progress logging, as well as some basic
  mathematical functions such as the rising factorial.
- `verify_and_refine_root.jl` - Implements the functions related to
  the interval Newton method. See the individual functions
  documentation for more details.
- `U.jl`, `U_expansion.jl` - Implements evaluation and expansion of
  the confluent hypergeometric function $U$. In particular this makes
  use of the bounds from Lemmas 7.1 and 7.2 in the self-similar blowup
  paper.
- `CGLParams.jl` - Defines the type `CGLParams` that is used for
  storing the parameters that are held fixed throughout the
  computations, that is $d$, $\omega$, $\sigma$ and $\delta$. It also
  contains code to generate the initial guesses for the NLS equation.
- `refine_approximation.jl` - Contains the code for non-rigorous
  refinement of initial approximations of zeros for the function $G$.
- `G_solve.jl` - Functions for computing roots of $G$ with the help of
  functions from `verify_and_refine_root.jl`.
- `G_hat_solve.jl` - Functions for computing roots of $\hat{G}$ with
  the help of functions from `verify_and_refine_root.jl`.
- `eigenvalues_finite_difference.jl` - Implementation of a finite
  difference method for the operator $L_{\hat{Q}}$ to compute initial
  approximations of eigenvalues.
- `src/CGLBranch.jl` - Contains the
  [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl)
  code used for computing numerical approximations of the branches and
  is completely self-contained (it is implemented as a submodule). It
  supports computing the branches using both $\epsilon$ and $\kappa$
  as continuation parameters, as well as two different approaches for
  approximating $Q_\infty$.
