#
# export JULIA_DEPOT_PATH=`pwd`/.julia
#
using Pkg
# Pkg.add( Pkg.PackageSpec(url="https://github.com/clima/climatemachine.jl", rev="master") )
Pkg.add( Pkg.PackageSpec(url="https://github.com/clima/climatemachine.jl", rev="78fd04b202a1c0c85c61bdc6bdf6c7f71218ec5c") )
Pkg.add("MPI")
Pkg.add("StaticArrays")
Pkg.activate("ClimateMachine")
Pkg.instantiate()

using ClimateMachine
