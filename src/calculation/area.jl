using Calculus
using Optim
using Formatting

const Φᵧ = 1.6180339887498948482

type DirtyNonadiabaticArea <: NonadiabaticArea
  states::Tuple{Int, Int}
  coordinate_from::Float64
  coordinate_to::Float64
  peaks::Array{Tuple{Float64, Float64}, 1}
  pits::Array{Tuple{Float64, Float64}, 1}
  sign::Int
  function DirtyNonadiabaticArea()
    this = new()
    this.states = (-1, -1)
    this.coordinate_from = -1
    this.coordinate_to = -1
    this.peaks = Array{Tuple{Float64, Float64}, 1}()
    this.pits = Array{Tuple{Float64, Float64}, 1}()
    this.sign = 0
    return this
  end
end

"""
NTRS (NASA Technical Reports Server)
Iott, J.; Haftka, R. T.; Adelman, H. M.
"Selecting step sizes in sensitivity analysis by finite differences"
https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19850025225.pdf
"""
function Δhₒₗ(ϵₐ, df_dx, d²f_dx²)
  Φ₀ = abs(d²f_dx²) ⋅ (1 + df_dx*df_dx)
  Φ = Φ₀ > 1e-4 ? Φ₀ : 1e-4
  return 2⋅√(ϵₐ / Φ)
end

function Δhₒₚₜ(ΔRₛₜ, df_dx)
  af = abs(df_dx)
  Δh₀ = ΔRₛₜ / (Φᵧ * (af < 1 ? 1 : af))
  return Δh₀ < ΔRₛₜ ? Δh₀ : ΔRₛₜ
end

function detectSinglePeakAreas(M_∂_∂R::Array{Function, 2}, nonadiabatic_config::NonadiabaticAreasConfiguration, Rstop::Float64)
  area_config = nonadiabatic_config.nonadiabatic_areas[SINGLE_PEAK::NonadiabaticAreaTypes]

  N = size(M_∂_∂R, 1)
  areas = Array{Array{NonadiabaticArea, 1}, 2}(N, N)
  fill!(areas, Array{NonadiabaticArea, 1}())

  Rₛₜₐᵣₜ = nonadiabatic_config.coordinate_start; ΔRₘₐₓ = nonadiabatic_config.coordinate_step; Rₛₜₒₚ = Rstop
  ΔRₚᵢₑₛₑ = nonadiabatic_config.coordinate_piece
  ϵₐ_y = abs(nonadiabatic_config.coordinate_step_error)

  ϵₚₑₐₖ = abs(area_config.error_∂_∂R_peak)
  yₛₘₐₗₗ = abs(area_config.vanishing_∂_∂R_value)
  ϵ_yₛₘₐₗₗ = abs(area_config.error_vanishing_∂_∂R_value)

  for i = 1:N, j = 1:N
    if i < j && j - i == 1
      dirty_areas = Array{DirtyNonadiabaticArea, 1}()
      areas[i, j] = Array{NonadiabaticArea, 1}()

      table = Array{Tuple{Float64, Float64}, 1}()
      τ = M_∂_∂R[i, j]
      R = Rₛₜₐᵣₜ
      while R <= Rₛₜₒₚ
        τᵥ = τ(R)
        dτ_dR = derivative(τ, R)
        push!(table, (R, τᵥ))
        ΔR = Δhₒₚₜ(ΔRₘₐₓ, dτ_dR)
        R += ΔR
      end

      # ***********
      # K. Villaverde, V. Kreinovich
      # "A linear-time algorithm that locates local extrema of
      # a function of one variable from interval measurement results"
      # https://www.researchgate.net/publication/244415993_A_linear-time_algorithm_that_locates_local_extrema_of_a_function_of_one_variable_from_interval_measurement_results
      # ***********
      n = size(table, 1)
      s = 0; M = table[1][2]; m = table[1][2]
      k = 1; xₖ = table[k][1]; yₖ = table[k][2]
      # -----------
      σₖ = sign(yₖ); χₖ = abs(abs(yₖ) - abs(yₛₘₐₗₗ))
      ℷ = false; ℵ = χₖ > ϵ_yₛₘₐₗₗ
      Α = nothing
      # -----------
      for k = 2:n
        xₖ = table[k][1]; yₖ = table[k][2]
        # -----------
        χₖ = abs(abs(yₖ) - abs(yₛₘₐₗₗ))
        if ℵ
          ℵ = χₖ > ϵ_yₛₘₐₗₗ
        else
          σₖ = sign(yₖ)
          Αₛₜₐᵣₜ = (χₖ > ϵ_yₛₘₐₗₗ && σₖ > 0 && s == 1) || (χₖ > ϵ_yₛₘₐₗₗ && σₖ < 0 && s == -1)
          Αₛₜₒₚ = (χₖ <= ϵ_yₛₘₐₗₗ && σₖ > 0 && s == -1) || (χₖ <= ϵ_yₛₘₐₗₗ && σₖ < 0 && s == 1)
          if Αₛₜₐᵣₜ
            if !ℷ
              ℷ = true
              Α = DirtyNonadiabaticArea()
              Α.states = (i, j)
              Α.coordinate_from = xₖ
              Α.sign = σₖ
            end
          elseif Αₛₜₒₚ
            if ℷ
              ℷ = false
              Α.coordinate_to = xₖ
              push!(dirty_areas, Α)
            end
          else
            # nothing
          end
        end
        # -----------
        if s == 0
          Mₑ = M - 2*ϵₚₑₐₖ; mₑ = m + 2*ϵₚₑₐₖ
          if Mₑ <= yₖ <= mₑ
            s = 0
          elseif Mₑ > yₖ
            s = -1
          elseif yₖ > mₑ
            s = 1
          end
          M = max(M, yₖ); m = min(m, yₖ)
        elseif s == 1
          Mₑ = M - 2*ϵₚₑₐₖ;
          if Mₑ <= yₖ
            M = max(M, yₖ)
            s = 1
          else
            for l in convert(Array{Int}, linspace(k-1, 1, k-1))
              xₗ = table[l][1]; yₗ = table[l][2]
              if yₗ < Mₑ
                # -----------
                if Α ≠ nothing
                  τᵢₙᵥ = R -> -τ(R)
                  result = Optim.optimize(τᵢₙᵥ, xₗ, xₖ, Optim.Brent())
                  xₘₐₓ = Optim.minimizer(result)
                  if Α.sign > 0
                    push!(Α.peaks, (xₘₐₓ, τ(xₘₐₓ)))
                  else
                    push!(Α.pits, (xₘₐₓ, τ(xₘₐₓ)))
                  end
                end
                # -----------
                #maximum = (xₗ, xₖ, M - ϵₚₑₐₖ)
                #push!(maxima, maximum)
                s = -1; m = yₖ
                break
              end
            end
          end
        elseif s == -1
          mₑ = m + 2*ϵₚₑₐₖ
          if mₑ >= yₖ
            m = min(m, yₖ)
            s = -1
          else
            for l in convert(Array{Int}, linspace(k-1, 1, k-1))
              xₗ = table[l][1]; yₗ = table[l][2]
              if yₗ > mₑ
                # -----------
                if Α ≠ nothing
                  result = Optim.optimize(τ, xₗ, xₖ, Optim.Brent())
                  xₘᵢₙ = Optim.minimizer(result)
                  if Α.sign > 0
                    push!(Α.pits, (xₘᵢₙ, τ(xₘᵢₙ)))
                  else
                    push!(Α.peaks, (xₘᵢₙ, τ(xₘᵢₙ)))
                  end
                end
                # -----------
                #minimum = (xₗ, xₖ, m - ϵₚₑₐₖ)
                #push!(minima, minimum)
                yₖ₁ = table[k-1][2]
                s = 1; M = yₖ₁
                break
              end
            end
          end
        else
          # empty
        end
      end

      if size(dirty_areas, 1) > 0
        single_peak_areas = Array{DirtyNonadiabaticArea, 1}()
        for darea in dirty_areas
          if size(darea.peaks, 1) == 1 && size(darea.pits, 1) == 0
            dpeak = first(darea.peaks)

            new_area = SinglePeakNonadiabaticArea()
            new_area.states = darea.states
            new_area.coordinate_∂_∂R = dpeak[1]
            new_area.value_∂_∂R = dpeak[2]
            new_area.coordinate_from = darea.coordinate_from
            new_area.coordinate_to = darea.coordinate_to
            new_area.sign = darea.sign
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
  M_Αˡᶻ = Array{Array{SinglePeakNonadiabaticArea, 1}, 2}(N, N)
  fill!(M_Αˡᶻ, Array{SinglePeakNonadiabaticArea, 1}())

  Α_config = nonadiabatic_config.nonadiabatic_areas[SINGLE_PEAK::NonadiabaticAreaTypes]
  ϵˡᶻ = abs(Α_config.error_potential_∂_∂R_coordinate)
  ϵₘᵢₙᵢₘᵤₘ = abs(Α_config.error_potential_distance_coordinate)

  for i = 1:N, j = 1:N
    if i >= j || abs(i - j) ≠ 1
      continue
    end
    M_Αˡᶻ[i, j] = Array{SinglePeakNonadiabaticArea, 1}()
    M_Αˡᶻ[j, i] = Array{SinglePeakNonadiabaticArea, 1}()

    Αₛ = areas[i, j]
    for Α in Αₛ
      Rₐ = Α.coordinate_from; Rᵦ = Α.coordinate_to
      R₀ = Α.coordinate_∂_∂R; τ = Α.value_∂_∂R; σ = Α.sign

      U₁(R) = M_Hₐ[i, i](R); U₂(R) = M_Hₐ[j, j](R);
      ΔU(R) = abs(U₂(R) - U₁(R))

      result = Optim.optimize(ΔU, Rₐ, Rᵦ, Optim.Brent())
      Rₘᵢₙᵢₘᵤₘ = Optim.minimizer(result)

      if abs(Rₘᵢₙᵢₘᵤₘ - R₀) <= ϵˡᶻ
        Αˡᶻ = SinglePeakNonadiabaticArea()
        Αˡᶻ.states = (i, j)
        Αˡᶻ.coordinate_∂_∂R = R₀; Αˡᶻ.value_∂_∂R = τ
        Αˡᶻ.coordinate_potentials = Rₘᵢₙᵢₘᵤₘ
        Αˡᶻ.coordinate_from = Rₐ; Αˡᶻ.coordinate_to = Rᵦ
        Αˡᶻ.sign = σ
        push!(M_Αˡᶻ[i, j], Αˡᶻ)

        Αˡᶻ_inv = SinglePeakNonadiabaticArea()
        Αˡᶻ_inv.states = (j, i)
        Αˡᶻ_inv.coordinate_∂_∂R = R₀; Αˡᶻ_inv.value_∂_∂R = -τ
        Αˡᶻ_inv.coordinate_potentials = Rₘᵢₙᵢₘᵤₘ
        Αˡᶻ_inv.coordinate_from = Rₐ; Αˡᶻ_inv.coordinate_to = Rᵦ
        Αˡᶻ_inv.sign = -σ
        push!(M_Αˡᶻ[j, i], Αˡᶻ_inv)
      else
        # nothing
      end
    end
  end

  return M_Αˡᶻ
end
