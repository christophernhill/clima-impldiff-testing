#
# export JULIA_DEPOT_PATH=`pwd`/.julia
#
using Pkg
Pkg.add( Pkg.PackageSpec(url="https://github.com/clima/climatemachine.jl", rev="master") )
Pkg.add("MPI")
Pkg.add("StaticArrays")
Pkg.activate("ClimateMachine")
Pkg.instantiate()

using ClimateMachine
