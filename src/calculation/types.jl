# -----------------------------------------------------------------------------
# Common Types
# -----------------------------------------------------------------------------
using DataFrames
using Formatting
import Base.show
import Base.-

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
function show(io::IO, Aᵛ::Vector{NonadiabaticArea})
  content = join(collect("$(A)" for A in Aᵛ), ",\n")
  print(io, "Α{$content}")
end
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
function -(A::SinglePeakNonadiabaticArea)
  Aᶜᵒⁿʲ = SinglePeakNonadiabaticArea()
  Aᶜᵒⁿʲ.states = reverse(A.states)
  Aᶜᵒⁿʲ.coordinate_∂_∂R = A.coordinate_∂_∂R
  Aᶜᵒⁿʲ.value_∂_∂R = -A.value_∂_∂R
  Aᶜᵒⁿʲ.coordinate_potentials = A.coordinate_potentials
  Aᶜᵒⁿʲ.coordinate_from = A.coordinate_from
  Aᶜᵒⁿʲ.coordinate_to = A.coordinate_to
  Aᶜᵒⁿʲ.sign = -A.sign
  return Aᶜᵒⁿʲ
end
function show(io::IO, A::SinglePeakNonadiabaticArea)
  i = A.states[1]; iˢᵘᵇ = int2indexsub(i); iᵐᵒˡ = int2molstate(i)
  j = A.states[2]; jˢᵘᵇ = int2indexsub(j); jᵐᵒˡ = int2molstate(j)
  sign = A.sign > 0 ? '+' : '-'
  str_R₀ᵁ = A.coordinate_potentials > 1e-10 ? format("{:.5f}", A.coordinate_potentials) : "undefined"
  const fe = FormatExpr("Α{1}{2}{3}[{4}{5} {6} {7}{8}|{9}{10}]{11}{12}={13:.5f}, {14}={15:.5f}, {16}={17:.5f}, {18}={19:.6e}, {20}={21}{22}")
  print(io,
    format(fe,
      iˢᵘᵇ, '⋅', jˢᵘᵇ, sign, '△', '→', '⟨', iᵐᵒˡ, jᵐᵒˡ, '⟩',
      '{',
        "Rᵃ", A.coordinate_from, "Rᵇ", A.coordinate_to,
        "R₀", A.coordinate_∂_∂R, "τ(R₀)", A.value_∂_∂R, "R₀ᵁ", str_R₀ᵁ,
      '}'))
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
function show(io::IO, A::LandauZenerArea)
  i = A.states[1]; iˢᵘᵇ = int2indexsub(i); iᵐᵒˡ = int2molstate(i)
  j = A.states[2]; jˢᵘᵇ = int2indexsub(j); jᵐᵒˡ = int2molstate(j)
  sign = A.τ₀ > 0 ? '+' : '-'
  const fe = FormatExpr("Α{1}{2}{3}[{4}{5} {6} {7}{8}|{9}{10}]{11}{12}={13:.5f}, {14}={15:.5f}, {16}={17:.5f}, {18}={19:.6e}, {20}:{21}{22}")
  print(io,
    format(fe,
      iˢᵘᵇ, '⋅', jˢᵘᵇ, sign, "Δˡᶻ", '→', '⟨', iᵐᵒˡ, jᵐᵒˡ, '⟩',
      '{',
        "Rₐ", A.Rₐ, "Rᵦ", A.Rᵦ, "R₀", A.R₀, "τ₀", A.τ₀, "∂/∂R", A.∂_∂R,
      '}'))
end
# ----------- non-adiabatic areas -----------

# ----------- calculation results -----------
# ------------------------------------------------------------------
