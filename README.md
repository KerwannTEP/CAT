# CAT
Chandrasekhar Anisotropic Theory

## Aim of the code

Computation of the NR diffusion coefficients used in the FP equation for an anisotropic Plummer model.

## Installation

Install Julia by following the instruction at `https://julialang.org/downloads/platform/`.

To invoke Julia in the Terminal, you need to make sure that the julia command-line program is in your `PATH`. 

On MacOS, we must create a link in `/usr/local/bin` (here for Julia 1.5):

```
$ sudo ln -s /Applications/Julia-1.5.app/Contents/Resources/julia/bin/julia /usr/local/bin/julia
```

## Julia packages (TODO)

Open the terminal in the folder `packages` and type

```
$ julia Install-pkg.jl
```

to install the following packages:

- `HDF5`
- `ArgParse`
- `HypergeometricFunctions`
- `SpecialFunctions`
- `StaticArrays`
- `Interpolations`
- `PolynomialRoots`

(TODO)

### !! WARNING !!

**DO NOT INTERRUPT THE DOWNLOADING OF THE PACKAGES !!!!**

The root finding algorithm employed in this library is described in

* J. Skowron & A. Gould, 2012, "General Complex Polynomial Root Solver and Its
  Further Optimization for Binary Microlenses",
  [arXiv:1203.1034](http://arxiv.org/abs/1203.1034)

This algorithm aims to be fast and precise, more than the well known `zroots`
procedure described in *Numerical Recipes* books, whose implementations in C and
Fortran are not available as free software, according to the
[definition](https://www.gnu.org/philosophy/free-sw.html) of the Free Software
Foundation.



## Compute the *LOCAL* diffusion coefficients in energy-angular momentum space

To compute the *LOCAL* NR diffusion coefficients in energy-momentum, one needs to open 
`code/compute/ComputeCoeffsEnergy.jl`.

Then, one needs to write for which radius `r`, energy `E`, angular momentum `L` and field star mass `m` one wants 
to compute the diffusion coefficients. Once this is done, save the file and run the command 

```
$ julia ComputeCoeffsEnergy.jl --q 0.0
```

where `--q 0.0` sets the cluster's anisotropy to `q=0.0`. This parameter needs to be a `Float64`.



## Compute the diffusion coefficients in action space

To compute the NR diffusion coefficients in action space `(Jr,L)`, one needs to open 
`code/compute/ComputeCoeffsAction.jl`.

Then, one needs to write for which radial action `Jr`, angular momentum `L` and field star mass `m`  one wants 
to compute the diffusion coefficients. Once this is done, save the file and run the command 

```
$ julia ComputeCoeffsAction.jl --q 0.0
```

where `--q 0.0` sets the cluster's anisotropy to `q=0.0`. This parameter needs to be a `Float64`.
