using MPI
using ClimateMachine
using ClimateMachine.Mesh.Topologies
using ClimateMachine.Mesh.Grids
using ClimateMachine.DGMethods
using ClimateMachine.DGMethods.NumericalFluxes

import ClimateMachine.DGMethods:
    init_state_auxiliary!, update_auxiliary_state!, update_auxiliary_state_gradient!, vars_state, VerticalDirection, init_state_prognostic!, boundary_state!, flux_first_order!, flux_second_order!, source!, wavespeed

import  ClimateMachine.SystemSolvers:
    BatchedGeneralizedMinimalResidual, linearsolve!

include("IVDCModel.jl")

 #####
 # Basic initialization
 #####
 ClimateMachine.init()
 ArrayType = ClimateMachine.array_type()
 mpicomm = MPI.COMM_WORLD
 FT = Float64

 #####
 # Create mesh
 #####
 const N = 4
 const Nˣ = 20
 const Nʸ = 20
 const Nᶻ = 20
 const Lˣ = 4e6  # m
 const Lʸ = 4e6  # m
 const H = 1000  # m

 xrange = range(FT(0); length = Nˣ + 1, stop = Lˣ)
 yrange = range(FT(0); length = Nʸ + 1, stop = Lʸ)
 zrange = range(FT(-H); length = Nᶻ + 1, stop = 0)

 brickrange_2D = (xrange, yrange)
 topl_2D       =
   BrickTopology( mpicomm, brickrange_2D, periodicity = (false, false) );
 grid_2D = DiscontinuousSpectralElementGrid(
   topl_2D,
   FloatType = FT,
   DeviceArray = ArrayType,
   polynomialorder = N,
 )
 brickrange_3D = (xrange, yrange, zrange)
 topl_3D = StackedBrickTopology(
   mpicomm,
   brickrange_3D;
   periodicity = (false, false, false),
   boundary = ((1, 1), (1, 1), (2, 3)),
 )
 grid_3D = DiscontinuousSpectralElementGrid(
   topl_3D,
   FloatType = FT,
   DeviceArray = ArrayType,
   polynomialorder = N,
 )

 # Wave speeds for numerical flux
 const cʰ = 1  # typical of ocean internal-wave speed
 const cᶻ = 0

 # Set a timestep for implicit solve
 dt = 5400

 # Create balance law and RHS arrays for diffusion equation
 ivdc_dg = IVDCDGModel(
  IVDCModel{FT}(;dt=3600,cʰ=cʰ,cᶻ=cᶻ),
  grid_3D,
  RusanovNumericalFlux(),
  CentralNumericalFluxSecondOrder(),
  CentralNumericalFluxGradient();
  direction=VerticalDirection(),
 );

 ivdc_Q   = init_ode_state(ivdc_dg, FT(0); init_on_cpu=true)
 ivdc_RHS = init_ode_state(ivdc_dg, FT(0); init_on_cpu=true)

 # Instantiate a batched GM-RES solver uaing balance law as its operator
 ivdc_bgm_solver=BatchedGeneralizedMinimalResidual(
  ivdc_dg,
  ivdc_Q;
  max_subspace_size=10);

 # Set up right hand side
 ivdc_Q.θ .= 0
 ivdc_RHS.θ .= 0
 ivdc_dg(ivdc_Q,ivdc_RHS,nothing,0;increment=false);

#     dg = OceanDGModel(
#        model,
#        grid_3D,
#    #   CentralNumericalFluxFirstOrder(),
#        RusanovNumericalFlux(),
#        CentralNumericalFluxSecondOrder(),
#        CentralNumericalFluxGradient(),
#    )




