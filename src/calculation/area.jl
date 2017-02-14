using Calculus
using Optim

type DirtyNonadiabaticArea <: NonadiabaticArea
  states::Tuple{Int, Int}
  coordinate_from::Float64
  coordinate_to::Float64
  peaks::Array{Tuple{Float64, Float64}, 1}
  pits::Array{Tuple{Float64, Float64}, 1}
  sign::Int
  function DirtyNonadiabaticArea()
    this = new()
    this.states = (-1 , -1)
    this.coordinate_from = -1
    this.coordinate_to = -1
    this.peaks = Array{Tuple{Float64, Float64}, 1}()
    this.pits = Array{Tuple{Float64, Float64}, 1}()
    this.sign = 0
  end
end

"""
NTRS (NASA Technical Reports Server)
Iott, J.; Haftka, R. T.; Adelman, H. M.
"Selecting step sizes in sensitivity analysis by finite differences"
https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19850025225.pdf
"""
function Δhₒₚₜ(ϵₐ, df_dx, d²f_dx²)
  Φ₀ = abs(d²f_dx²) ⋅ (1 + df_dx*df_dx)
  Φ = Φ₀ > 1e-8 ? Φ₀ : 1e-8
  return 2⋅√(ϵₐ / Φ)
end

function detectSinglePeakAreas(M_∂_∂R::Array{Function, 2}, nonadiabatic_config::NonadiabaticAreasConfiguration, Rstop::Float64)
  area_config = nonadiabatic_config.nonadiabatic_areas[SINGLE_PEAK::NonadiabaticAreaTypes]

  N = size(M_∂_∂R, 1)
  areas = Array{Array{NonadiabaticArea, 1}, 2}(N, N)
  fill!(areas, Array{NonadiabaticArea, 1}())

  Rₛₜₐᵣₜ = nonadiabatic_config.coordinate_start; ΔRₘₐₓ = nonadiabatic_config.coordinate_step; Rₛₜₒₚ = Rstop
  ϵₐ_τ = abs(nonadiabatic_config.coordinate_step_error)

  ϵₚₑₐₖ = abs(area_config.error_∂_∂R_peak)
  τₛₘₐₗₗ = abs(area_config.vanishing_∂_∂R_value)
  ϵ_τₛₘₐₗₗ =abs(area_config.error_vanishing_∂_∂R_value)

  for i = 1:N, j = 1:N
    if i < j && j - i == 1
      dirty_areas = Array{DirtyNonadiabaticArea, 1}()

      # ------------------------------------
      # ----------- Before Cycle -----------
      τ = M_∂_∂R[i, j]
      τᵥ = τ(Rₛₜₐᵣₜ); dτ_dR = derivative(τ, Rₛₜₐᵣₜ); d²τ_dR² = second_derivative(τ, Rₛₜₐᵣₜ)

      ΔRₒₚₜ = Δhₒₚₜ(ϵₐ_τ, dτ_dR, d²τ_dR²)
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
          ΔRₒₚₜ = Δhₒₚₜ(ϵₐ_τ, dτ_dR, d²τ_dR²)
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

      current_area = DirtyNonadiabaticArea()
      current_area.states = (i, j)
      # ----------- Cycle -----------
      ℵ = false
      R = Rₛₜₐᵣₜ+ΔR
      while R ≤ Rₛₜₒₚ
        τᵥ_ΔR = τ(R); dτ_dR_ΔR = derivative(τ, R); d²τ_dR²_ΔR = second_derivative(τ, R)
        # -----------

        τₕₑᵢᵧₕₜ = abs(abs(τᵥ) - τₛₘₐₗₗ)
        if τₕₑᵢᵧₕₜ < ϵ_τₛₘₐₗₗ
          if ℵ
            current_area.coordinate_to = R
            if size(current_area, 1) == 0
              throw(ErrorException("Incorrect area definition for <$i|∂/∂R|$j> at [$(current_area.coordinate_from), $(current_area.coordinate_to)]."))
            end
            some_peak = first(current_area.peaks)
            current_area.sign = sign(some_peak[2])

            push!(dirty_areas, current_area)

            current_area = DirtyNonadiabaticArea()

            ℵ = false
          end
        elseif τₕₑᵢᵧₕₜ ≥ ϵ_τₛₘₐₗₗ
          if !ℵ
            current_area.coordinate_from = R
          end
          ℵ = true

          sgn_τᵥ = sign(τᵥ); sgn_τᵥ_ΔR = sign(τᵥ_ΔR)
          sgn_dτ_dR = sign(dτ_dR); sgn_dτ_dR_ΔR = sign(dτ_dR_ΔR)
          sgn_d²τ_dR² = sign(d²τ_dR²); sgn_d²τ_dR²_ΔR = sign(d²τ_dR²_ΔR)

          Δτ = abs(τᵥ_ΔR - τᵥ)
          if sgn_dτ_dR == -sgn_dτ_dR_ΔR
            if sgn_τᵥ == -sgn_τᵥ_ΔR
              throw(ErrorException("Unexpected situation for <$i|∂/∂R|$j> at [$(R - ΔR), $R]; τᵥ=$τᵥ, τᵥ_ΔR=$τᵥ_ΔR."))
            end
            #if sgn_d²τ_dR² == -sgn_d²τ_dR²_ΔR
            #  throw(ErrorException("Unexpected situation for second derivative of <$i|∂/∂R|$j> at [$(R - ΔR), $R]; d²τ_dR²=$d²τ_dR², d²τ_dR²_ΔR_ΔR=$d²τ_dR²_ΔR."))
            #end

            Rₐ = R - ΔR; Rᵦ = R
            Rₑₓₜᵣₑₘᵤₘ = R; τₑₓₜᵣₑₘᵤₘ = τᵥ
            if sgn_d²τ_dR² < 0 || (sgn_dτ_dR > 0 && sgn_dτ_dR_ΔR < 0)
              itfunc(x) = -τ(x)
              result = Optim.optimize(itfunc, Rₐ, Rᵦ, Optim.Brent())
              if !Optim.converged(result) || Optim.abs_tol > ϵₚₑₐₖ
                throw(ErrorException("Unable to find a local maximum for <$i|∂/∂R|$j> at [$(R - ΔR), $R]; τᵥ=$τᵥ, τᵥ_ΔR=$τᵥ_ΔR; converged=$(Optim.converged(result)); ϵₐᵦₛ=$(Optim.abs_tol(result))."))
              end
              Rₑₓₜᵣₑₘᵤₘ = Optim.minimum(result); τₑₓₜᵣₑₘᵤₘ = τ(Rₑₓₜᵣₑₘᵤₘ)

              if sgn_τᵥ > 0
                push!(current_area.peaks, (Rₑₓₜᵣₑₘᵤₘ, τₑₓₜᵣₑₘᵤₘ))
              else
                push!(current_area.pits, (Rₑₓₜᵣₑₘᵤₘ, τₑₓₜᵣₑₘᵤₘ))
              end
            elseif sgn_d²τ_dR² > 0 || (sgn_dτ_dR < 0 && sgn_dτ_dR_ΔR > 0)
              tfunc(x) = τ(x)
              result = Optim.optimize(tfunc, Rₐ, Rᵦ, Optim.Brent())
              if !Optim.converged(result) || Optim.abs_tol > ϵₚₑₐₖ
                throw(ErrorException("Unable to find a local minimum for <$i|∂/∂R|$j> at [$(R - ΔR), $R]; τᵥ=$τᵥ, τᵥ_ΔR=$τᵥ_ΔR; converged=$(Optim.converged(result)); ϵₐᵦₛ=$(Optim.abs_tol(result))."))
              end
              Rₑₓₜᵣₑₘᵤₘ = Optim.minimum(result); τₑₓₜᵣₑₘᵤₘ = τ(Rₑₓₜᵣₑₘᵤₘ)

              if sgn_τᵥ > 0
                push!(current_area.pits, (Rₑₓₜᵣₑₘᵤₘ, τₑₓₜᵣₑₘᵤₘ))
              else
                push!(current_area.peaks, (Rₑₓₜᵣₑₘᵤₘ, τₑₓₜᵣₑₘᵤₘ))
              end
            else
              # do nothing
            end
          end
        else
          throw(ErrorException("Unexpected extremum search execution point for <$i|∂/∂R|$j> at [$(R - ΔR), $R]; τᵥ=$τᵥ, τᵥ_ΔR=$τᵥ_ΔR."))
        end

        # -----------
        τᵥ = τᵥ_ΔR; dτ_dR = dτ_dR_ΔR; d²τ_dR² = d²τ_dR²_ΔR
        ΔRₒₚₜ = Δhₒₚₜ(ϵₐ_τ, dτ_dR, d²τ_dR²)
        ΔR = ΔRₘₐₓ > ΔRₒₚₜ ? ΔRₒₚₜ : ΔRₘₐₓ
        # -----------
        R += ΔR
      end
      # ----------- End of Cycle -----------
      # ------------------------------------

      if size(dirty_areas, 1) > 0
        single_peak_areas = Array{DirtyNonadiabaticArea, 1}()
        for darea in dirty_areas
          if size(darea.peaks, 1) == 1 && size(darea.pits, 1) == 0
            dpeak = first(darea.peaks)

            new_area = SinglePeakNonadiabaticArea()
            new_area.states = darea.states
            new_area.coordinate_∂_∂R = dpeak[1]
            new_area.coordinate_from = darea.coordinate_from
            new_area.coordinate_to = darea.coordinate_to
            push!(areas[i, j], new_area)
          end
        end
      end
    else
      # undef reference
    end
  end

  for i = 1:N, j = 1:N
    if i > j && i - j == 1
      areas[i, j] = areas[j, i]
    end
  end

  return areas
end

function detectLandauZenerAreas(M_Hₐ::Array{Function, 2}, areas::Array{Array{NonadiabaticArea, 1}, 2}, nonadiabatic_config::NonadiabaticAreasConfiguration, Rstop::Float64)
  N = size(M_Hₐ, 1)
  areas_lz = Array{Array{SinglePeakNonadiabaticArea, 1}, 2}(N, N)

  area_config = nonadiabatic_config.nonadiabatic_areas[SINGLE_PEAK::NonadiabaticAreaTypes]
  ϵ_lz = abs(area_config.error_potential_∂_∂R_coordinate)
  ϵₘᵢₙᵢₘᵤₘ = abs(area_config.error_potential_distance_coordinate)

  for i = 1:N, j = 1:N
    if i >= j && abs(i - j) ≠ 1
      continue
    end

    candidate = areas[i, j]
    Rₐ = candidate.coordinate_from; Rᵦ = candidate.coordinate_to
    R₀ = candidate.coordinate_∂_∂R

    U₁(R) = M_Hₐ[i, i](R); U₂(R) = M_Hₐ[j, j](R);
    ΔU(R) = abs(U₂(R) - U₁(R))

    result = Optim.optimize(ΔU, Rₐ, Rᵦ, Optim.Brent())
    Rₘᵢₙᵢₘᵤₘ = Optim.minimum(result)

    if abs(Rₘᵢₙᵢₘᵤₘ - R₀) <= ϵ_lz
      area_lz = SinglePeakNonadiabaticArea()
      area_lz.states = (i, j)
      area_lz.coordinate_∂_∂R = R₀; area_lz.coordinate_potentials = Rₘᵢₙᵢₘᵤₘ
      area_lz.coordinate_from = Rₐ; area_lz.coordinate_to = Rᵦ
      areas_lz[i, j] = area_lz

      area_lz_inv = SinglePeakNonadiabaticArea()
      area_lz_inv.states = (j, i)
      area_lz_inv.coordinate_∂_∂R = R₀; area_lz_inv.coordinate_potentials = Rₘᵢₙᵢₘᵤₘ
      area_lz_inv.coordinate_from = Rₐ; area_lz_inv.coordinate_to = Rᵦ
      areas_lz[j, i] = area_lz_inv
    end

    return areas_lz
  end
end
