# Standalone solver test with CLiMA driver

## One time setup in clean directory, dedicated project tree and with local MPI install in /usr/local
```
export JULIA_DEPOT_PATH=`pwd`/.julia
/Applications/Julia-1.4.app/Contents/Resources/julia/bin/julia --project=ClimateMachine mysetup.jl
```

## Run
```
export JULIA_DEPOT_PATH=`pwd`/.julia
/Applications/Julia-1.4.app/Contents/Resources/julia/bin/julia --project=ClimateMachine driver.jl
```

