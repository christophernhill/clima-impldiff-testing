# Linear model equations, for ocean model implicit vertical diffusion
# convective adjustment step.
#
using ClimateMachine.DGMethods.NumericalFluxes:
    NumericalFluxFirstOrder, NumericalFluxSecondOrder, NumericalFluxGradient

using ClimateMachine.DGMethods.NumericalFluxes:
    CentralNumericalFluxGradient, CentralNumericalFluxSecondOrder

newstyle=:(
using ClimateMachine.BalanceLaws:
   BalanceLaw, Prognostic, Auxiliary, Gradient, GradientFlux
)
oldstyle=:(
using ClimateMachine.BalanceLaws:
    BalanceLaw
)
@oldnew

using ClimateMachine.VariableTemplates

newstyle=:(
import ClimateMachine.DGMethods:
     init_state_auxiliary!, update_auxiliary_state!, update_auxiliary_state_gradient!, vars_state, VerticalDirection, boundary_state!
)
oldstyle=:(
 import ClimateMachine.DGMethods:
    init_state_auxiliary!, update_auxiliary_state!, update_auxiliary_state_gradient!, VerticalDirection, boundary_state!;
 import ClimateMachine.DGMethods:
    vars_state_auxiliary, vars_state_gradient, vars_state_conservative, vars_state_gradient_flux
)
@oldnew

using StaticArrays

"""
 IVDCModel{M} <: BalanceLaw

 This code defines DG `BalanceLaw` terms for an operator, L, that is evaluated from iterative
 implicit solver to solve an equation of the form

 (L + 1/Δt) ϕ^{n+1} = ϕ^{n}/Δt

  where L is a vertical diffusion operator with a spatially varying diffusion
  coefficient.

 # Usage

"""

# Create a new linear model instance
abstract type AbstractIVDCModel <: BalanceLaw end
struct IVDCModel{FT} <: AbstractIVDCModel 
 dt::FT
 κᶻ::FT
 κᶜ::FT
 cʰ::FT
 cᶻ::FT
 function IVDCModel{FT}(
  ;
  dt = FT(5400),
  κᶻ = FT(1e-4),
  κᶜ = FT(0.1),
  cʰ = FT(0),
  cᶻ = FT(0),
 ) where {FT<: AbstractFloat}
    return new{FT}(
      dt,
      κᶻ,
      κᶜ,
      cʰ,
      cᶻ,
    )
 end
end

function IVDCDGModel(
   bl::IVDCModel,
   grid,
   nfnondiff,
   nfdiff,
   gnf;
   kwargs...,
   )

   modeldata=(dt=bl.dt,κᶻ=bl.κᶻ,κᶜ=bl.κᶜ,cʰ=bl.cʰ,cᶻ=bl.cᶻ)

   return DGModel(bl,grid,nfnondiff,nfdiff,gnf;kwargs...,modeldata=modeldata,)
end

"""
 Set model state variables and operators
"""

# State variable and initial value, just one for now, θ
##
oldstyle=:(
vars_state_conservative(m::IVDCModel, FT) = @vars(θ::FT)
)
newstyle=:(
vars_state(m::IVDCModel, ::Prognostic, FT) = @vars(θ::FT)
)
@oldnew

newstyle=:(
init_state_prognostic!( args...) = ( init_state_cp!( args... ) )
)
oldstyle=:(
init_state_conservative!( args...) = ( init_state_cp!( args... ) )
)
@oldnew
function init_state_cp!( m::IVDCModel, Q::Vars, A::Vars, coords, t,)
  @inbounds begin
    Q.θ = -0
  end
  return nothing
end
##

##
oldstyle=:(
vars_state_auxiliary(m::IVDCModel, FT) = @vars(θ_init::FT)
)
newstyle=:(
vars_state(m::IVDCModel, ::Auxiliary, FT) = @vars(θ_init::FT)
)
@oldnew

function init_state_auxiliary!(m::IVDCModel,A::Vars, _...)
  @inbounds begin
    A.θ_init = -0
  end
  return nothing
end
##
#

# Variables and operations used in differentiating first derivatives
##
oldstyle=:(
vars_state_gradient(m::IVDCModel, FT) = @vars(∇θ::FT, ∇θ_init::FT,)
)
newstyle=:(
vars_state(m::IVDCModel, ::Gradient, FT) = @vars(∇θ::FT, ∇θ_init::FT,)
)
@oldnew

@inline function compute_gradient_argument!(
    m::IVDCModel,
    G::Vars,
    Q::Vars,
    A,
    t,
)
    #NAN G.∇θ = Q.θ
    #NAN G.∇θ_init = A.θ_init

    return nothing
end
##
#

# Variables and operations used in differentiating second derivatives
##
oldstyle=:(
vars_state_gradient_flux(m::IVDCModel, FT) = @vars(κ∇θ::SVector{3, FT})
)
newstyle=:(
vars_state(m::IVDCModel, ::GradientFlux, FT) = @vars(κ∇θ::SVector{3, FT})
)
@oldnew

@inline function compute_gradient_flux!(
    m::IVDCModel,
    D::Vars,
    G::Grad,
    Q::Vars,
    A::Vars,
    t,
)
    κ = diffusivity_tensor(m, G.∇θ_init[3])
    D.κ∇θ = -κ * G.∇θ .* false
    return nothing
end
##

##
## Set vertical diffusivity profile based on vertical hydrography profile
@inline function diffusivity_tensor(m::IVDCModel, ∂θ∂z)
    κᶻ = m.κᶻ
    κᶜ = m.κᶜ
    ∂θ∂z < 0 ? κ = (@SVector [0, 0, κᶜ]) : κ = (@SVector [0, 0, κᶻ])
    κ = (@SVector [0, 0, 0])
    return Diagonal(-κ.*0)
end
##
#

# Function to apply I to state variable

##
@inline function source!(
    m::IVDCModel,
    S::Vars,
    Q::Vars,
    D::Vars,
    A::Vars,
    t,
    direction,
)
    ivdc_dt = m.dt
    @inbounds begin
     S.θ = Q.θ/ivdc_dt
    end

    return nothing
end
##
#

# Numerical fluxes and boundaries

##
function flux_first_order!(::IVDCModel, _...) end

function flux_second_order!(
    ::IVDCModel,
    F::Grad,
    S::Vars,
    D::Vars,
    H::Vars,
    A::Vars,
    t,
)
    #NAN F.θ += D.κ∇θ
    F.θ = D.κ∇θ
end

function wavespeed(m::IVDCModel, n⁻, _...)
    # C = abs(SVector(m.cʰ, m.cʰ, m.cᶻ)' * n⁻)
    C = abs(SVector(0, 0, 0)' * n⁻)
    return C
end
##

function boundary_state!(
    nf::Union{NumericalFluxFirstOrder, NumericalFluxGradient, CentralNumericalFluxGradient},
    m::IVDCModel,
    Q⁺,
    A⁺,
    n,
    Q⁻,
    A⁻,
    bctype,
    t,
    _...,
)
    Q⁺.θ = Q⁻.θ
    return nothing
end

###    From -  function numerical_boundary_flux_gradient! , DGMethods/NumericalFluxes.jl
###    boundary_state!(
###        numerical_flux,
###        balance_law,
###        state_conservative⁺,
###        state_auxiliary⁺,
###        normal_vector,
###        state_conservative⁻,
###        state_auxiliary⁻,
###        bctype,
###        t,
###        state1⁻,
###        aux1⁻,
###    )

function boundary_state!(
    nf::Union{NumericalFluxSecondOrder,CentralNumericalFluxSecondOrder},
    m::IVDCModel,
    Q⁺,
    D⁺,
    A⁺,
    n⁻,
    Q⁻,
    D⁻,
    A⁻,
    bctype,
    t,
    _...,
)
    Q⁺.θ = Q⁻.θ
    D⁺.κ∇θ = n⁻ * -0 .* false
    return nothing
end

###    boundary_state!(
###        numerical_flux,
###        balance_law,
###        state_conservative⁺,
###        state_gradient_flux⁺,
###        state_auxiliary⁺,
###        normal_vector,
###        state_conservative⁻,
###        state_gradient_flux⁻,
###        state_auxiliary⁻,
###        bctype,
###        t,
###        state1⁻,
###        diff1⁻,
###        aux1⁻,
###    )

###    boundary_flux_second_order!(
###        numerical_flux,
###        balance_law,
###        Grad{S}(flux),
###        state_conservative⁺,
###        state_gradient_flux⁺,
###        state_hyperdiffusive⁺,
###        state_auxiliary⁺,
###        normal_vector,
###        state_conservative⁻,
###        state_gradient_flux⁻,
###        state_hyperdiffusive⁻,
###        state_auxiliary⁻,
###        bctype,
###        t,
###        state1⁻,
###        diff1⁻,
###        aux1⁻,
###    )

