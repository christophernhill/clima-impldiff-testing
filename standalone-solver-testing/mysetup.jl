#
# export JULIA_DEPOT_PATH=`pwd`/.julia
#
using Pkg
Pkg.add( Pkg.PackageSpec(url="https://github.com/clima/climatemachine.jl", rev="master") )
ENV["JULIA_MPI_DIR"]="/usr/local"
Pkg.add("MPI")
Pkg.add("StaticArrays")
Pkg.activate("ClimateMachine")
Pkg.instantiate()

using ClimateMachine
