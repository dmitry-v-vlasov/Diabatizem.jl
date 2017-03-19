using Calculus
using Optim
using Formatting
using Logging

type DirtyNonadiabaticArea <: NonadiabaticArea
  states::Tuple{Int, Int}
  coordinate_from::Float64
  coordinate_to::Float64
  peaks::Vector{Tuple{Float64, Float64}}
  pits::Vector{Tuple{Float64, Float64}}
  sign::Int
  function DirtyNonadiabaticArea()
    this = new()
    this.states = (-1, -1)
    this.coordinate_from = -1
    this.coordinate_to = -1
    this.peaks = Vector{Tuple{Float64, Float64}}()
    this.pits = Vector{Tuple{Float64, Float64}}()
    this.sign = 0
    return this
  end
end

function detectSinglePeakAreas(M_∂_∂R::Array{Function, 2}, M_∂_∂Rᵈᵃᵗᵃ::Array{Float64, 2}, nonadiabatic_config::NonadiabaticAreasConfiguration, Rstop::Float64)
  Logging.configure(level=INFO)

  Rₛₜₐᵣₜ = nonadiabatic_config.coordinate_start; ΔRₘₐₓ = nonadiabatic_config.coordinate_step; Rₛₜₒₚ = Rstop
  ΔRₚᵢₑₛₑ = nonadiabatic_config.coordinate_piece
  ϵₐ_y = abs(nonadiabatic_config.coordinate_step_error)

  # -----------
  M_∂_∂R_sorted = sortrows(M_∂_∂Rᵈᵃᵗᵃ; by=row->(row[1]))
  M_∂_∂R_vector_filtered = Vector{Vector{Float64}}()
  for row in IteratorRow(M_∂_∂R_sorted)
    R = row[1]
    if Rₛₜₐᵣₜ <= R <= Rₛₜₒₚ
      push!(M_∂_∂R_vector_filtered, row)
    end
  end
  L = size(M_∂_∂R_vector_filtered, 1); Nᶜ = size(M_∂_∂R_vector_filtered[1], 1)
  M_∂_∂Rᵍᵒᵒᵈ = Array{Float64, 2}(L, Nᶜ)
  for l = 1:L
    M_∂_∂Rᵍᵒᵒᵈ[l, :] = M_∂_∂R_vector_filtered[l]
  end
  # -----------

  area_config = nonadiabatic_config.nonadiabatic_areas[SINGLE_PEAK::NonadiabaticAreaTypes]
  ϵₚₑₐₖ = abs(area_config.error_∂_∂R_peak)
  yₛₘₐₗₗ = abs(area_config.vanishing_∂_∂R_value)
  ϵ_yₛₘₐₗₗ = abs(area_config.error_vanishing_∂_∂R_value)

  N = size(M_∂_∂R, 1)
  areas = Array{Vector{NonadiabaticArea}, 2}(N, N)
  fill!(areas, Vector{NonadiabaticArea}())

  info("Single peak non-adiabatic area detection.")
  info("Search configuration: [$Rₛₜₐᵣₜ, $Rₛₜₒₚ], ΔRₘₐₓ=$ΔRₘₐₓ; R in data table: [$(M_∂_∂Rᵍᵒᵒᵈ[1, 1]), $(M_∂_∂Rᵍᵒᵒᵈ[L, 1])]; ϵ(⟨|∂/∂R|⟩ₚₑₐₖ)=$ϵₚₑₐₖ; ⟨|∂/∂R|⟩ₛₘₐₗₗ=$yₛₘₐₗₗ; ϵ(⟨|∂/∂R|⟩ₛₘₐₗₗ)=$ϵ_yₛₘₐₗₗ")
  for i = 1:N, j = 1:N
    if i < j && j - i == 1
      info("⇩⇩⇩⇩⇩⇩⇩⇩⇩⇩⇩ Scanning ⟨$i|∂/∂R|$j⟩ ⇩⇩⇩⇩⇩⇩⇩⇩⇩⇩⇩")
      τ = M_∂_∂R[i, j]

      dirty_areas = Vector{DirtyNonadiabaticArea}()
      areas[i, j] = Vector{NonadiabaticArea}()

      table = Vector{Tuple{Float64, Float64}}()
      for row in IteratorRow(M_∂_∂Rᵍᵒᵒᵈ)
        R = row[1]; τᵗ = row[dataColumnOfSymetricMatrix(i, j, N) + 1]
        push!(table, (R, τᵗ))
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
      σₖ = sign(yₖ); χₖ = abs(yₖ) > abs(yₛₘₐₗₗ) ? abs(abs(yₖ) - abs(yₛₘₐₗₗ)) : abs(yₖ)
      ℷ = false; ℵ = χₖ > ϵ_yₛₘₐₗₗ
      Α = nothing
      # -----------
      for k = 2:n
        xₖ = table[k][1]; yₖ = table[k][2]
        # -----------
        χₖ = abs(yₖ) > abs(yₛₘₐₗₗ) ? abs(abs(yₖ) - abs(yₛₘₐₗₗ)) : abs(yₖ)
        if ℵ
          ℵ = χₖ > ϵ_yₛₘₐₗₗ
          #info("ℵ was true and now ℵ=$ℵ; χₖ=$χₖ; ϵ_yₛₘₐₗₗ=$ϵ_yₛₘₐₗₗ; xₖ=$xₖ; yₖ=$yₖ")
        else
          σₖ = sign(yₖ)
          Αₛₜₐᵣₜ = (χₖ > ϵ_yₛₘₐₗₗ && σₖ > 0 && s == 1) || (χₖ > ϵ_yₛₘₐₗₗ && σₖ < 0 && s == -1)
          Αₛₜₒₚ = (χₖ <= ϵ_yₛₘₐₗₗ && σₖ > 0 && s == -1) || (χₖ <= ϵ_yₛₘₐₗₗ && σₖ < 0 && s == 1)
          #info("ℵ was false anf now ℵ=$ℵ; χₖ=$χₖ; ϵ_yₛₘₐₗₗ=$ϵ_yₛₘₐₗₗ; σₖ=$σₖ; s=$s; ℷ=$ℷ; Αₛₜₐᵣₜ=$Αₛₜₐᵣₜ; Αₛₜₒₚ=$Αₛₜₒₚ; xₖ=$xₖ; yₖ=$yₖ")
          if Αₛₜₐᵣₜ
            if !ℷ
              ℷ = true
              Α = DirtyNonadiabaticArea()
              Α.states = (i, j)
              Α.coordinate_from = xₖ
              Α.sign = σₖ
              info("Αₛₜₐᵣₜ at $(Α.coordinate_from); χₖ=$χₖ; yₖ=$yₖ; ℵ=$ℵ; ℷ=$ℷ")
            end
          elseif Αₛₜₒₚ
            if ℷ
              ℷ = false
              Α.coordinate_to = xₖ
              push!(dirty_areas, Α)
              info("Αₛₜₒₚ at $(Α.coordinate_to); χₖ=$χₖ; yₖ=$yₖ; ℵ=$ℵ; ℷ=$ℷ")
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
            for l in convert(Vector{Int}, linspace(k-1, 1, k-1))
              xₗ = table[l][1]; yₗ = table[l][2]
              if yₗ < Mₑ
                # -----------
                if ℷ
                  τᵢₙᵥ = R -> -τ(R)
                  result = Optim.optimize(τᵢₙᵥ, xₗ, xₖ, Optim.Brent())
                  xₘₐₓ = Optim.minimizer(result)
                  τₘₐₓ = -Optim.minimum(result)
                  if Α.sign > 0
                    push!(Α.peaks, (xₘₐₓ, τₘₐₓ))
                  else
                    push!(Α.pits, (xₘₐₓ, τₘₐₓ))
                  end
                  info("𝔐aximum of ⟨$i|∂/∂R|$j⟩ found: Rₘₐₓ=$(format("{:.5f}", xₘₐₓ)); τ(Rₘₐₓ)=$(format("{:.6e}", τₘₐₓ)); ϵʳᵉˡ=$(format("{:.6e}", Optim.rel_tol(result))); ϵᵃᵇˢ=$(format("{:.6e}", Optim.abs_tol(result)))")
                end
                # -----------
                #info("ALG: Found function maximum in an area; τ(xₘₐₓ)≈$(M - ϵₚₑₐₖ); xₗ=$xₗ; xₖ=$xₖ")
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
            for l in convert(Vector{Int}, linspace(k-1, 1, k-1))
              xₗ = table[l][1]; yₗ = table[l][2]
              if yₗ > mₑ
                # -----------
                if ℷ
                  result = Optim.optimize(τ, xₗ, xₖ, Optim.Brent())
                  xₘᵢₙ = Optim.minimizer(result)
                  τₘᵢₙ = Optim.minimum(result)
                  if Α.sign > 0
                    push!(Α.pits, (xₘᵢₙ, τₘᵢₙ))
                  else
                    push!(Α.peaks, (xₘᵢₙ, τₘᵢₙ))
                  end
                  info("𝔐inimum of ⟨$i|∂/∂R|$j⟩ found: Rₘᵢₙ=$(format("{:.5f}", xₘᵢₙ)); τ(Rₘᵢₙ)=$(format("{:.6e}", τₘᵢₙ)); ϵʳᵉˡ=$(format("{:.6e}", Optim.rel_tol(result))); ϵᵃᵇˢ=$(format("{:.6e}", Optim.abs_tol(result)))")
                end
                # -----------
                #info("ALG: Found function minimum in an area; τ(xₘᵢₙ)≈$(M - ϵₚₑₐₖ); xₗ=$xₗ; xₖ=$xₖ")
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
        single_peak_areas = Vector{DirtyNonadiabaticArea}()
        for darea in dirty_areas
          if size(darea.peaks, 1) == 1 && size(darea.pits, 1) == 0
            dpeak = first(darea.peaks)
            R₀ = dpeak[1]; τₚₑₐₖ = dpeak[2]
            χₚₑₐₖ = abs(τₚₑₐₖ) > abs(yₛₘₐₗₗ) ? abs(abs(τₚₑₐₖ) - abs(yₛₘₐₗₗ)) : abs(τₚₑₐₖ)

            new_area = SinglePeakNonadiabaticArea()
            new_area.states = darea.states
            new_area.coordinate_∂_∂R = dpeak[1]
            new_area.value_∂_∂R = dpeak[2]
            new_area.coordinate_from = darea.coordinate_from
            new_area.coordinate_to = darea.coordinate_to
            new_area.sign = darea.sign
            new_area.coordinate_potentials = 0.0
            if χₚₑₐₖ <= ϵ_yₛₘₐₗₗ
              warn("Skipping low peak (χₚₑₐₖ=$χₚₑₐₖ ≤ ϵ(⟨|∂/∂R|⟩ₛₘₐₗₗ)=$ϵ_yₛₘₐₗₗ): $(new_area)")
              continue
            end
            push!(areas[i, j], new_area)
          end
        end
      end
      info("⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧⇧")
    else
      # undef reference
    end
  end

  for i = 1:N, j = 1:N
    if i > j && i - j == 1
      if (!isdefined(areas[j, i]) || isempty(areas[j, i])) continue end
      K = length(areas[j, i])
      conjugate_areas = Vector{NonadiabaticArea}(0)
      setindex!(areas, i, j, conjugate_areas)
      for k = 1:K
        push!(conjugate_areas, -(areas[j, i][k]))
      end
      areas[i, j] = conjugate_areas
    end
  end

  info()
  info("▫▫▫▫▫▫▫▫▫▫▫ Single Peak Areas Summary ▫▫▫▫▫▫▫▫▫▫▫")
  for i = 1:N, j=1:N
    if i < j && j - i == 1
      info("*********** Areas of ⟨$i|∂/∂R|$j⟩ ***********")
      for Aᵛ in areas[i, j]
        info(Aᵛ)
      end
      info("*******************************************")
    end
  end
  info("▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫▫")

  return areas
end

function detectLandauZenerAreas(M_Hₐ::Array{Function, 2}, areas::Array{Vector{NonadiabaticArea}, 2}, nonadiabatic_config::NonadiabaticAreasConfiguration, Rstop::Float64)
  N = size(M_Hₐ, 1)
  M_Αˡᶻ = Array{Vector{SinglePeakNonadiabaticArea}, 2}(N, N)
  fill!(M_Αˡᶻ, Vector{SinglePeakNonadiabaticArea}())

  Α_config = nonadiabatic_config.nonadiabatic_areas[SINGLE_PEAK::NonadiabaticAreaTypes]
  ϵˡᶻ = abs(Α_config.error_potential_∂_∂R_coordinate)
  ϵₘᵢₙᵢₘᵤₘ = abs(Α_config.error_potential_distance_coordinate)

  for i = 1:N, j = 1:N
    if i >= j || abs(i - j) ≠ 1
      continue
    end
    M_Αˡᶻ[i, j] = Vector{SinglePeakNonadiabaticArea}()
    M_Αˡᶻ[j, i] = Vector{SinglePeakNonadiabaticArea}()

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
