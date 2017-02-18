# -----------------------------------------------------------------------------
# Common Types
# -----------------------------------------------------------------------------

using DataFrames

# ------------------------------------------------------------------
type ParameterMatrix
  distance::Float64
  data::Array{Float64, Float64}
  function ParameterMatrix() new() end
end
# ------------------------------------------------------------------

# ----------- calculation results -----------
abstract Calculation

# ----------- asymptotics -----------
type PotentialAsymptotic <: Calculation
  state::Int
  value::Float64
  coordinate::Float64
  PotentialAsymptotic() = new(-1, -1, -1)
  function PotentialAsymptotic(a_state::Int, a_value::Float64, a_coordinate::Float64)
    this = new()
    this.state = a_state
    this.value = a_value
    this.coordinate = a_coordinate
    return this
  end
end

type AsymptoticCalculation <: Calculation
  potentials::Array{PotentialAsymptotic}
end
# ----------- asymptotics -----------

# ----------- non-adiabatic areas -----------
abstract NonadiabaticArea <: Calculation
type SinglePeakNonadiabaticArea <: NonadiabaticArea
  states::Tuple{Int, Int}
  coordinate_∂_∂R::Float64
  value_∂_∂R::Float64
  coordinate_potentials::Float64
  coordinate_from::Float64
  coordinate_to::Float64
  sign::Int
  function SinglePeakNonadiabaticArea() new() end
end

type LandauZenerArea <: NonadiabaticArea
  states::Tuple{Int, Int}
  R₀::Float64
  τ₀::Float64
  Rₐ::Float64
  Rᵦ::Float64
  ∂_∂R::Function
  function LandauZenerArea(instates::Tuple{Int, Int}, inR₀::Float64, inτ₀::Float64, inRₐ::Float64, inRᵦ::Float64)
    this = new()
    this.states = instates
    this.R₀ = inR₀; this.τ₀ = inτ₀
    this.Rₐ = inRₐ; this.Rᵦ = inRᵦ
    this.∂_∂R = R -> this.τ₀ / ((R - this.R₀)^2 + 4*(this.τ₀)^2)
    return this
  end
end
# ----------- non-adiabatic areas -----------

# ----------- diabatic data -----------
type DiabaticData
  hamiltonian::Dict{Float64, ParameterMatrix}
  ∂_∂R::Dict{Float64, ParameterMatrix}
  function DiabaticData() new() end
end

type DiabaticDataOutput
  hamiltonian::DataFrame
  ∂_∂R::DataFrame
  function DiabaticDataOutput() new() end
end
# ----------- diabatic data -----------

# ----------- all calculation results -----------
type GeneralCalculation <: Calculation
  asymptotics::AsymptoticCalculation
  nonadiabatic_areas::Array{NonadiabaticArea}
  function GeneralCalculation() new() end
end
# ----------- all calculation results -----------

# ----------- calculation results -----------
# ------------------------------------------------------------------
