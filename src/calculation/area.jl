using Calculus

type DirtyNonadiabaticArea <: NonadiabaticArea
  states::Tuple{Int, Int}
  coordinate_from::Float64
  coordinate_to::Float64
  peaks::Array{Tuple{Float64, Float64}, 1}
  pits::Array{Tuple{Float64, Float64}, 1}
  function DirtyNonadiabaticArea()
    this = new()
    this.states = (-1 , -1)
    this.peaks = Array{Tuple{Float64, Float64}, 1}()
    this.pits = Array{Tuple{Float64, Float64}, 1}()
  end
end

function detectSinglePeakAreas(M_∂_∂R::Array{Function, 2}, nonadiabatic_config::NonadiabaticAreasConfiguration, Rstop::Float64)
  N = size(M_∂_∂R, 1)
  areas = Array{SinglePeakNonadiabaticArea, 1}(N, N)

  R_start = nonadiabatic_areas.coordinate_start; ΔR = nonadiabatic_areas.coordinate_step; R_stop = Rstop

  area_config = nonadiabatic_config.nonadiabatic_areas[SINGLE_PEAK::NonadiabaticAreaTypes]

  ϵ_peak = area_config.error_∂_∂R_peak
  τ_small = area_config.vanishing_∂_∂R_value
  ϵ_τ_small = area_config.error_vanishing_∂_∂R_value

  for i = 1:N, j = 1:N
    if i < j && j - i == 1
      points = Float64[]
      τ = M_∂_∂R[i, j]
      τᵥ = τ(R_start); ∂τ_∂R = derivative(τ, R_start); ∂²τ_∂R² = second_derivative(τ, R_start)
      for R = R_start+ΔR:ΔR:R_stop
        τᵥ_ΔR = τ(R); ∂τ_∂R_ΔR = derivative(τ, R); ∂²τ_∂R²_ΔR = second_derivative(τ, R)
        # -----------

        sgn_∂τ_∂R = sign(∂τ_∂R); sgn_∂τ_∂R_ΔR = sign(∂τ_∂R_ΔR)
        #sgn_∂²τ_∂R² = sign(∂²τ_∂R²); sgn_∂²τ_∂R²_ΔR = sign(∂²τ_∂R²_ΔR)

        Δτ = abs(τᵥ_ΔR - τᵥ)
        if sgn_∂τ_∂R == -sgn_∂τ_∂R_ΔR
          if Δτ < ϵ_peak
            push!(points, (R + (R + ΔR))/2)
            R += ΔR
          else
            R -= ΔR
            ΔR /= 2
          end
        else
          if Δτ < ϵ_peak
            ΔR *= 2
          end
        end

        # -----------
        τᵥ = τᵥ_ΔR; ∂τ_∂R = ∂τ_∂R_ΔR; ∂²τ_∂R² = ∂²τ_∂R²_ΔR
      end
    else
      # undef reference
    end
  end
end
