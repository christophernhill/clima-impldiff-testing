# Standalone solver test with CLiMA driver

## One time setup in clean directory, dedicated project tree and with local MPI install in /usr/local
## Assume ClimateMachine.jl is in clone one level up
```
export JULIA_DEPOT_PATH=`pwd`/.julia
export JULIA_MPI_DIR="/usr/local"
/Applications/Julia-1.4.app/Contents/Resources/julia/bin/julia --project=../ClimateMachine.jl mysetup.jl
```

## Run
```
export JULIA_DEPOT_PATH=`pwd`/.julia
export JULIA_MPI_DIR="/usr/local"
/Applications/Julia-1.4.app/Contents/Resources/julia/bin/julia --project=../ClimateMachine.jl driver.jl
```

