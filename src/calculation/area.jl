using Calculus
using Optim
using Formatting

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
  Δh₀ = ΔRₛₜ / (af < 1 ? 1 : af)
  return Δh₀ < ΔRₛₜ ? Δh₀ : ΔRₛₜ
end

function detectSinglePeakAreas(M_∂_∂R::Array{Function, 2}, nonadiabatic_config::NonadiabaticAreasConfiguration, Rstop::Float64)
  area_config = nonadiabatic_config.nonadiabatic_areas[SINGLE_PEAK::NonadiabaticAreaTypes]

  N = size(M_∂_∂R, 1)
  areas = Array{Array{NonadiabaticArea, 1}, 2}(N, N)
  fill!(areas, Array{NonadiabaticArea, 1}())

  Rₛₜₐᵣₜ = nonadiabatic_config.coordinate_start; ΔRₘₐₓ = nonadiabatic_config.coordinate_step; Rₛₜₒₚ = Rstop
  ΔRₚᵢₑₛₑ = nonadiabatic_config.coordinate_piece
  ϵₐ_τ = abs(nonadiabatic_config.coordinate_step_error)

  ϵₚₑₐₖ = abs(area_config.error_∂_∂R_peak)
  τₛₘₐₗₗ = abs(area_config.vanishing_∂_∂R_value)
  ϵ_τₛₘₐₗₗ = abs(area_config.error_vanishing_∂_∂R_value)

  for i = 1:N, j = 1:N
    if i < j && j - i == 1
      dirty_areas = Array{DirtyNonadiabaticArea, 1}()

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

      maxima = Array{Tuple{Float64, Float64, Float64}, 1}()
      minima = Array{Tuple{Float64, Float64, Float64}, 1}()

      n = size(table, 1)
      s = 0; M = table[1][2]; m = table[1][2]
      xᵢ = table[i][1]; yᵢ = table[i][2]
      for i = 2:n
        xᵢ = table[i][1]; yᵢ = table[i][2]
        if s == 0
          Mₑ = M - 2*ϵₚₑₐₖ; mₑ = m + 2*ϵₚₑₐₖ
          if Mₑ <= yᵢ <= mₑ
            s = 0
          elseif Mₑ > yᵢ
            s = -1
          elseif yᵢ > mₑ
            s = 1
          end
          M = max(M, yᵢ); m = min(m, yᵢ)
        elseif s == 1
          Mₑ = M - 2*ϵₚₑₐₖ;
          if Mₑ <= yᵢ
            M = max(M, yᵢ)
            s = 1
          else
            for j in convert(Array{Int}, linspace(i-1, 1, i-1))
              xⱼ = table[j][1]; yⱼ = table[j][2]
              if yⱼ < Mₑ
                maximum = (xⱼ, xᵢ, M - ϵₚₑₐₖ)
                push!(maxima, maximum)
                s = -1; m = yᵢ
                break
              end
            end
          end
        elseif s == -1
          mₑ = m + 2*ϵₚₑₐₖ
          if mₑ >= yᵢ
            m = min(m, yᵢ)
            s = -1
          else
            for j in convert(Array{Int}, linspace(i-1, 1, i-1))
              xⱼ = table[j][1]; yⱼ = table[j][2]
              if yⱼ > mₑ
                minimum = (xⱼ, xᵢ, m - ϵₚₑₐₖ)
                push!(minima, minimum)
                yᵢ₁ = table[i-1][2]
                s = 1; M = yᵢ₁
                break
              end
            end
          end
        else
          # empty
        end
      end

      println(maxima)
      println("#####")
      println(minima)
      return

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
