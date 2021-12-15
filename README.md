# CAT
Chandrasekhar Anisotropic Theory

## Aim of the code

Computation of the non-resonant (NR) diffusion coefficients used in the orbit-averaged Fokker-Planck equation for an anisotropic Plummer model.

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



## Compute the **local** diffusion coefficients in energy-angular momentum space

To compute the **local** NR diffusion coefficients in energy-momentum, one needs to open 
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



## Compute a diffusion coefficients map in action space

To compute an action space map of the NR diffusion coefficients for a cluster with `10^5` single-mass stars,
one needs to access the `code/compute` folder and run the following command in 
the terminal:

```
$ julia ActionOrbitalMap.jl --parallel no --q 0.0
```

where `--q 0.0` sets the cluster's anisotropy to `q=0.0` and `-- parallel no` runs the code without parallelisation.

If one wants to run this with parallelisation, one needs to run the following 
commands (supposing one is using bash)

```
$ export JULIA_NUM_THREADS=12
$ export JULIA_CPU_THREADS=12
$ julia -p 12 ActionOrbitalMap.jl --parallel yes --q 0.0
```
	
where 12 is the number of parallelised threads. One can check the number of 
threads by opening the Julia terminal and by running the command

```
julia> Threads.nthreads()
```

The resulting file will be created in the folder `code/data` under the name 
`Dump_Diffusion_Coefficients_Action_Orbital_Map_q_0.0.hf5`, where the suffix `_q_0.0.` depends
on the anisotropy parameter given in argument `--q`.

This `.hf5` files contains the following tables:

- `q` : value of the anisotropy parameter of the cluster.
- `nbJrMeasure`: number of sampling loci in radial action.
- `nbLMeasure`: number of sampling loci in angular momentum.
- `tabL` : 1D-table of angular momentum values used for the evaluation of the coefficients.
- `tabJr` : 1D-table of radial action values used for the evaluation of the coefficients.
- `tabLJr`: 1D-table of the action space loci `(L,Jr)` used for the evaluation of the coefficients.
- `tabDNRJr`: 1D-table of the coefficients DNR_Jr corresponding to the action space loci `(L,Jr)` in table `tabLJr`.
- `tabDNRJr`: 1D-table of the coefficients DNR_L corresponding to the action space loci `(L,Jr)` in table `tabLJr`.
- `tabDNRJrJr`: 1D-table of the coefficients DNR_JrJr corresponding to the action space loci `(L,Jr)` in table `tabLJr`.
- `tabDNRJrL`: 1D-table of the coefficients DNR_JrL corresponding to the action space loci `(L,Jr)` in table `tabLJr`.
- `tabDNRLL`: 1D-table of the coefficients DNR_LL corresponding to the action space loci `(L,Jr)` in table `tabLJr`.


Note that the sampling in action space is done in log-log. One can modify the region of sampling by differently modifying the section `Action space parameter` in the file `ActionOrbitalMap.jl`.

Those files can be recovered using, for example, `Mathematica`, through the command

```
Import[NotebookDirectory[] <> StringJoin["data/Dump_Diffusion_Coefficients_Action_Orbital_Map_q_0.0.hf5"], {"Datasets", "tabDNRJr"}]
```
here for the (already computed )table `tabDNRJr` for a cluster with anisotropy `q=0.0`.




## Compute the diffusion coefficients in action space

To compute the NR diffusion coefficients in action space `(Jr,L)`, one needs to open 
`code/compute/ComputeCoeffsAction.jl`.

Then, one needs to write for which radial action `Jr`, angular momentum `L` and field star mass `m`  one wants 
to compute the diffusion coefficients. Once this is done, save the file and run the command 

```
$ julia ComputeCoeffsAction.jl --q 0.0
```

where `--q 0.0` sets the cluster's anisotropy to `q=0.0`. This parameter needs to be a `Float64`.



## Compute dF/dt in action space

To compute an action space map of `dF/dt` for a cluster with `10^5` single-mass stars,
one needs to access the `code/compute` folder and run the following command in 
the terminal:

```
$ julia MappingdFdt.jl --parallel no --q 0.0
```

where `--q 0.0` sets the cluster's anisotropy to `q=0.0` and `-- parallel no` runs the code without parallelisation.

If one wants to run this with parallelisation, one needs to run the following 
commands (supposing one is using bash)

```
$ export JULIA_NUM_THREADS=12
$ export JULIA_CPU_THREADS=12
$ julia -p 12 MappingdFdt.jl --parallel yes --q 0.0
```
	
where 12 is the number of parallelised threads. One can check the number of 
threads by opening the Julia terminal and by running the command

```
julia> Threads.nthreads()
```

The resulting file will be created in the folder `code/data` under the name 
`Dump_dFdt_Map_q_0.0.hf5`, where the suffix `_q_0.0.` depends
on the anisotropy parameter given in argument `--q`.

This `.hf5` files contains the following tables:

- `q` : value of the anisotropy parameter of the cluster.
- `nbJrMeasure`: number of sampling loci in radial action.
- `nbLMeasure`: number of sampling loci in angular momentum.
- `Mtot`: the total mass of the cluster.
- `Npart`: the number of stars in the cluster.
- `tabL` : 1D-table of angular momentum values used for the evaluation of the coefficients.
- `tabJr` : 1D-table of radial action values used for the evaluation of the coefficients.
- `tabLJr`: 1D-table of the action space loci `(L,Jr)` used for the evaluation of the coefficients.
- `tabdFdt`: 1D-table of the values of dF/dt corresponding to the action space loci `(L,Jr)` in table `tabLJr`.


Note that the sampling in action space is done linearly. One can modify the region of sampling by differently modifying the section `Action space parameter` in the file `ActionOrbitalMap.jl`.

Those files can be recovered using, for example, `Mathematica`, through the command

```
Import[NotebookDirectory[] <> StringJoin["data/Dump_dFdt_Map_q_0.0.hf5"], {"Datasets", "tabdFdt"}]
```
here for the (already computed )table `tabdFdt` for a cluster with anisotropy `q=0.0`.

