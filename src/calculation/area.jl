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
    this.coordinate_from = -1
    this.coordinate_to = -1
    this.peaks = Array{Tuple{Float64, Float64}, 1}()
    this.pits = Array{Tuple{Float64, Float64}, 1}()
  end
end

function detectSinglePeakAreas(M_∂_∂R::Array{Function, 2}, nonadiabatic_config::NonadiabaticAreasConfiguration, Rstop::Float64)
  area_config = nonadiabatic_config.nonadiabatic_areas[SINGLE_PEAK::NonadiabaticAreaTypes]

  N = size(M_∂_∂R, 1)
  areas = Array{SinglePeakNonadiabaticArea, 1}(N, N)

  Rₛₜₐᵣₜ = nonadiabatic_config.coordinate_start; ΔRₘₐₓ = nonadiabatic_config.coordinate_step; Rₛₜₒₚ = Rstop
  ϵₐ_τ = abs(nonadiabatic_config.coordinate_step_error)

  ϵₚₑₐₖ = abs(area_config.error_∂_∂R_peak)
  τₛₘₐₗₗ = abs(area_config.vanishing_∂_∂R_value)
  ϵ_τₛₘₐₗₗ =abs(area_config.error_vanishing_∂_∂R_value)

  for i = 1:N, j = 1:N
    if i < j && j - i == 1
      dirtyAreas = Array{DirtyNonadiabaticArea, 1}()

      # ------------------------------------
      # ----------- Before Cycle -----------
      τ = M_∂_∂R[i, j]
      τᵥ = τ(Rₛₜₐᵣₜ); dτ_dR = derivative(τ, Rₛₜₐᵣₜ); d²τ_dR² = second_derivative(τ, Rₛₜₐᵣₜ)

      # NTRS (NASA Technical Reports Server)
      # Iott, J.; Haftka, R. T.; Adelman, H. M.
      # "Selecting step sizes in sensitivity analysis by finite differences"
      # https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19850025225.pdf
      ΔRₒₚₜ = 2⋅√(ϵₐ_τ/(abs(d²τ_dR²) > 1e-8 ? abs(d²τ_dR²) : 1e-8))
      ΔR = ΔRₘₐₓ > ΔRₒₚₜ ? ΔRₒₚₜ : ΔRₘₐₓ

      # -----------
      # -----------
      if abs(abs(τᵥ) - τₛₘₐₗₗ) > ϵ_τₛₘₐₗₗ
        R = Rₛₜₐᵣₜ+ΔR
        while R ≤ Rₛₜₒₚ
          τᵥ_ΔR = τ(R); dτ_dR_ΔR = derivative(τ, R); d²τ_dR²_ΔR = second_derivative(τ, R)
          # -----------
          τₕₑᵢᵧₕₜ = abs(abs(τᵥ) - τₛₘₐₗₗ)
          if abs(abs(τᵥ) - τₛₘₐₗₗ) < ϵ_τₛₘₐₗₗ
            Rₛₜₐᵣₜ = R
            break
          end
          # -----------
          τᵥ = τᵥ_ΔR; dτ_dR = dτ_dR_ΔR; d²τ_dR² = d²τ_dR²_ΔR
          ΔRₒₚₜ = 2⋅√(ϵₐ_τ/(abs(d²τ_dR²) < 1e-8 ? abs(d²τ_dR²) : 1e-8))
          ΔR = ΔRₘₐₓ > ΔRₒₚₜ ? ΔRₒₚₜ : ΔRₘₐₓ
          # -----------
          R += ΔR
        end
        if Rₛₜₐᵣₜ >= Rₛₜₒₚ || abs(Rₛₜₒₚ - Rₛₜₐᵣₜ) < ϵₚₑₐₖ
          throw(ErrorException("Unable to find a zero point to start in <$i|∂/∂R|$j>."))
        end
      end
      # -----------
      # -----------

      currentArea = DirtyNonadiabaticArea()
      currentArea.states = (i, j)
      # ----------- Cycle -----------
      ℵ = false
      R = Rₛₜₐᵣₜ+ΔR
      while R ≤ Rₛₜₒₚ
        τᵥ_ΔR = τ(R); dτ_dR_ΔR = derivative(τ, R); d²τ_dR²_ΔR = second_derivative(τ, R)
        # -----------

        τₕₑᵢᵧₕₜ = abs(abs(τᵥ) - τₛₘₐₗₗ)
        if τₕₑᵢᵧₕₜ < ϵ_τₛₘₐₗₗ
          if ℵ
            currentArea.coordinate_to = R
            areas[i, j] = currentArea
            currentArea = DirtyNonadiabaticArea()

            ℵ = false
          end
        elseif τₕₑᵢᵧₕₜ ≥ ϵ_τₛₘₐₗₗ
          ℵ = true

          sgn_τᵥ = sign(τᵥ); sgn_τᵥ_ΔR = sign(τᵥ_ΔR)
          sgn_dτ_dR = sign(dτ_dR); sgn_dτ_dR_ΔR = sign(dτ_dR_ΔR)
          sgn_d²τ_dR² = sign(d²τ_dR²); sgn_d²τ_dR²_ΔR = sign(d²τ_dR²_ΔR)

          Δτ = abs(τᵥ_ΔR - τᵥ)
          if sgn_dτ_dR == -sgn_dτ_dR_ΔR
            if sgn_τᵥ == -sgn_τᵥ_ΔR
              throw(ErrorException("Unexpected situation for <$i|∂/∂R|$j> at [$(R - ΔR), $R]."))
            end
            if sgn_d²τ_dR² == -sgn_d²τ_dR²_ΔR
              throw(ErrorException("Unexpected situation for second derivative of <$i|∂/∂R|$j> at [$(R - ΔR), $R]."))
            end

            Rₑₓₜᵣₑₘᵤₘ = R; τₑₓₜᵣₑₘᵤₘ = τᵥ
            if sgn_d²τ_dR² < 0
              if sgn_τᵥ > 0
                push!(currentArea.peaks, (Rₑₓₜᵣₑₘᵤₘ, τₑₓₜᵣₑₘᵤₘ))
              else
                push!(currentArea.pits, (Rₑₓₜᵣₑₘᵤₘ, τₑₓₜᵣₑₘᵤₘ))
              end
            else
              if sgn_τᵥ > 0
                push!(currentArea.pits, (Rₑₓₜᵣₑₘᵤₘ, τₑₓₜᵣₑₘᵤₘ))
              else
                push!(currentArea.peaks, (Rₑₓₜᵣₑₘᵤₘ, τₑₓₜᵣₑₘᵤₘ))
              end
            end
          end
        else
          # nothing
        end

        # -----------
        τᵥ = τᵥ_ΔR; dτ_dR = dτ_dR_ΔR; d²τ_dR² = d²τ_dR²_ΔR
        ΔRₒₚₜ = 2⋅√(ϵₐ_τ/(abs(d²τ_dR²) < 1e-8 ? abs(d²τ_dR²) : 1e-8))
        ΔR = ΔRₘₐₓ > ΔRₒₚₜ ? ΔRₒₚₜ : ΔRₘₐₓ
        # -----------
        R += ΔR
      end
      # ----------- End of Cycle -----------
      # ------------------------------------
    else
      # undef reference
    end
  end
end
