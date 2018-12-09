using Calculus
using Logging
using ProgressMeter
using Combinatorics
using Nullables

using PolynomialRoots

import Dierckx
import NumericalMath

function diabatize(
    Hₐ::Array{Function, 2},
    ∂_∂R::Array{Function, 2}, ∂_∂Rᵐᵒᵈᵉˡ::Array{Function, 2}, ∂_∂R_arg::Vector{Float64},
    solutions::Vector{LocalSolution}, C::DiabatizationSettings, LZ::Array{Array{LandauZenerArea, 1}, 2})
    Logging.configure(level=INFO)
    info("**************************************************")
    info("Diabatizing...")
    R = unique(sort(collect(Base.flatten(map(solution->solution.points, solutions)))))
    info("Derived grid: [$(R[1]) ... $(R[end])] ($(length(R)))")
    @assert issorted(R)
    steps = R[2:end] - R[1:end-1]
    h = minimum(abs(steps)); H = maximum(abs(steps))
    info("Diabatizing in the interval [$(R[1]), $(R[end])] with hₘᵢₙ=$h, hₘₐₓ=$H")
    Rᶜ = clearGrid(R, 1e-10)
    empty!(steps)
    steps = Rᶜ[2:end] - Rᶜ[1:end-1]
    h = minimum(abs(steps)); H = maximum(abs(steps))
    info("Cleared grid [$(Rᶜ[1]), $(Rᶜ[end])] ($(length(Rᶜ))) with hₘᵢₙ=$h, hₘₐₓ=$H")
    info("Making partial matrices")
    Sᶠ = Vector{Array{Function, 2}}()
    Sᶠˢᵖ = Vector{Array{Dierckx.Spline1D, 2}}()
    for solution ∈ solutions
        Sᶠⁱ = matl2matfsl(solution.points, solution.S)
        push!(Sᶠ, Sᶠⁱ[1])
        push!(Sᶠˢᵖ, Sᶠⁱ[2])
    end
    info("Number of partial transformation matrices length - $(length(Sᶠ))")
    Rᶜ, Sᵛᵉᶜ, Hᵈ, ∂_∂Rᵈ, ∂_∂Rᵐ = diabatizeWithPartialMatrices(Hₐ, ∂_∂R, ∂_∂Rᵐᵒᵈᵉˡ, ∂_∂R_arg, Rᶜ, Sᶠ, Sᶠˢᵖ, solutions, LZ)
    info("**************************************************")
    return Rᶜ, Sᵛᵉᶜ, Hᵈ, ∂_∂Rᵈ, ∂_∂Rᵐ
end

function diabatizeWithPartialMatrices(
    Hₐ::Array{Function, 2},
    ∂_∂R::Array{Function, 2}, ∂_∂Rᵐᵒᵈᵉˡ::Array{Function, 2}, ∂_∂R_arg::Vector{Float64},
    Rᶜ::Vector{Float64}, Sᶠ::Vector{Array{Function, 2}}, Sᶠˢᵖ::Vector{Array{Dierckx.Spline1D, 2}},
    Sl::Vector{LocalSolution}, LZ::Array{Array{LandauZenerArea, 1}, 2})
    Nˡ = length(Sᶠ)
    N = size(Sᶠ[1], 1)
    @assert N == size(Hₐ, 1)
    # ----
    info("Full transformation matrix computation...")
    Sᵛ = Vector{Array{Float64, 2}}()
    for R ∈ Rᶜ
        S = eye(N, N)
        for i = 1:Nˡ
            Sⁱ = matf2mat(R, Sᶠ[i])
            S = S * Sⁱ
        end
        push!(Sᵛ, S)
    end
    info("Done.")
    # ----
    Sᵛᶠ, Sᵛᶠˢᵖ = matl2matfsl(Rᶜ, Sᵛ)
    # ----
    info("Transforming...")
    Hᵈ = Vector{Array{Float64, 2}}()
    ∂_∂Rᵈ = Vector{Array{Float64, 2}}()
    ∂_∂Rᵐ = Vector{Array{Float64, 2}}()
    Sᵛᵉᶜ = Vector{Array{Float64, 2}}()

    interval_states = Dict{Vector{Int}, Tuple{Float64, Float64}}()
    foreach(s -> begin interval_states[s.states] = s.interval end, Sl)
    info("Interval states: $interval_states")
    info("Making new potentials...")
    for R ∈ Rᶜ
        S = matf2mat(R, Sᵛᶠ)
        S⁻¹ = S'

        Hᴬ = matf2mat(R, Hₐ)
        Hᴰ = S⁻¹*Hᴬ*S

        @assert typeof(Hᴰ) == Array{Float64, 2}
        @assert size(Hᴰ, 1) == N

        push!(Hᵈ, Hᴰ)
        push!(Sᵛᵉᶜ, S)
    end

    info("Making new ⟨·|∂/∂R|·⟩...")
    for R ∈ ∂_∂R_arg
        S = matf2mat(R, Sᵛᶠ)
        ∇S = Dierckx.derivative.(Sᵛᶠˢᵖ, R; nu=1)
        S⁻¹ = S'

        ∂_∂Rᴬ = matf2mat(R, ∂_∂R); ∂_∂Rᴹ = matf2mat(R, ∂_∂Rᵐᵒᵈᵉˡ)
        ∂_∂Rᴰ = S⁻¹ * ∂_∂Rᴬ * S + S⁻¹ * ∇S

        state_list = map(is -> is[1], filter(is -> is[2][1] ≤ R ≤ is[2][2], collect(interval_states)))
        for states ∈ state_list
            @assert(issorted(states), "The states $states are not sorted.")
            s¹ = states[1]; sᵉ = states[end]
            #∂_∂Rᴰ[s¹:sᵉ, s¹:sᵉ] = S⁻¹[s¹:sᵉ, s¹:sᵉ] * (∂_∂Rᴬ[s¹:sᵉ, s¹:sᵉ] - ∂_∂Rᴹ[s¹:sᵉ, s¹:sᵉ]) * S[s¹:sᵉ, s¹:sᵉ]
            ∂_∂Rᴰ[s¹:sᵉ, s¹:sᵉ] = zeros(sᵉ-s¹+1, sᵉ-s¹+1)
            ∂_∂Rᴰ[s¹:sᵉ, s¹:sᵉ] = ∂_∂Rᴬ[s¹:sᵉ, s¹:sᵉ] - ∂_∂Rᴹ[s¹:sᵉ, s¹:sᵉ]
        end

        @assert typeof(∂_∂Rᴰ) == Array{Float64, 2}
        push!(∂_∂Rᵈ, ∂_∂Rᴰ)
        push!(∂_∂Rᵐ, ∂_∂Rᴹ)
    end

    # sols = Dict{Vector{Int}, LocalSolution}()
    # foreach(s -> begin sols[s.states] = s end, Sl)
    ∂_∂Rᴰᵈᵃᵗᵃ = matl2matlupperx(∂_∂Rᵈ);
    N_lz = size(LZ, 1)
    @assert size(LZ, 1) == size(LZ, 2)
    @assert N_lz == N
    info("=================================Smoothing ⟨i|∂/∂R|j⟩=====================================")
    for i = 1:N_lz, j = 1:N_lz
        if i >= j
            continue
        end
        if !isempty(LZ[i, j]) && i < j
            info("********** Smoothing areas for ⟨$i|∂/∂R|$j⟩ **********")
            for area ∈ LZ[i, j]
                for sol ∈ Sl
                    info("----------")
                    states = sol.states
                    s¹ = minimum(states); sᵉ = maximum(states)
                    states_ddr = combinations(collect(s¹:sᵉ), 2)
                    for states_ddr_ij ∈ states_ddr
                        R₀ = area.R₀
                        R¹ = area.Rₐ; Rᵉ = area.Rᵦ
                        ΔR = abs(Rᵉ - R¹)
                        R¹ = R¹ - ΔR; Rᵉ = Rᵉ + ΔR;
                        peaks = sol.peaks
                        peak_found_at_R₀ = findfirst(peak->abs(peak[1] - R₀) < 1e-1, peaks) > 0
                        if (i == states_ddr_ij[1] && j == states_ddr_ij[2]) && peak_found_at_R₀
                            info("!!!! - START smoothing ⟨$i|∂/∂R|$j⟩ - !!!!")
                            @assert all(pair->pair[1] < pair[2], states_ddr)
                            states_ddr_f = filter(states -> states[2] - states[1] == 1, combinations(collect(s¹:sᵉ), 2))

                            info("Smoothing area: $(LZ[i, j]) in interval [$R¹, $Rᵉ]")
                            l¹ = findlast(R -> R < R¹, ∂_∂R_arg); lᵉ = findlast(R -> R <= Rᵉ, ∂_∂R_arg)
                            info("Curve smoothing for the solution:\n$sol\nin interval [$R¹, $Rᵉ] for ⟨$i|∂/∂R|$j⟩...")
                            k = dataColumnOfSymetricMatrix(i, j, N)
                            vR = ∂_∂R_arg[l¹:lᵉ]
                            ddr_sample = ∂_∂Rᴰᵈᵃᵗᵃ[l¹:lᵉ, k]

                            It = NumericalMath.trapz(vR, ddr_sample)
                            α = 4 / (Rᵉ - R¹)
                            β = - 2 / tan(2 * It - π)
                            γ = - (Rᵉ - R₀) * (R₀ - R¹) / (Rᵉ - R¹)
                            @assert β^2 - 4 * α * γ > 0

                            τ₁₂ = real(roots([γ, β, α]))

                            info("Roots for [$R¹, - $R₀ - , $Rᵉ]: $(τ₁₂)")
                            abs_max = findmax(abs(τ₁₂))
                            τ₀ = τ₁₂[abs_max[2]]

                            ddr_spectrum = fft(ddr_sample)
                            Lˢ = length(ddr_spectrum)
                            Lᵉᵈᵍᵉ = floor(Int, Lˢ / 4)
                            ΔLˢ = abs(Lˢ - Lᵉᵈᵍᵉ)
                            high_f = filter(f -> abs(f) > 5.0, abs(ddr_spectrum[Lᵉᵈᵍᵉ:end]))
                            mean_f = mean(abs(ddr_spectrum[Lᵉᵈᵍᵉ:end]))
                            if length(high_f) / ΔLˢ > 0.4 && mean_f > 5.0
                                ϵ_τ = 0.005
                                ΔRˢᵐ = 2 * abs(τ₀) * √((1 - ϵ_τ) / ϵ_τ)
                                R¹ˢ = R₀ - ΔRˢᵐ; Rᵉˢ = R₀ + ΔRˢᵐ
                                info("Recalculated smoothing area: ⟨$i|∂/∂R|$j⟩ in interval [$R¹ˢ, $Rᵉˢ]")
                                l¹ˢ = findlast(R -> R < R¹ˢ, ∂_∂R_arg); lᵉˢ = findlast(R -> R <= Rᵉˢ, ∂_∂R_arg)
                                vR = ∂_∂R_arg[l¹ˢ:lᵉˢ]

                                info("For the coupling ⟨$i|∂/∂R|$j⟩ found R₀ = $R₀ and τ₀ = $τ₀, ∫⟨$i|∂/∂R|$j⟩dR = $It, [$R¹ˢ, $Rᵉˢ], ΔRˢᵐ=$ΔRˢᵐ")
                                τ(R) = τ₀ / ((R - R₀)^2 + 4 * τ₀ * τ₀)
                                ddr_sample_new = τ.(vR)
                                ∂_∂Rᴰᵈᵃᵗᵃ[l¹ˢ:lᵉˢ, k] = ddr_sample_new
                                info("!!!! - END smoothing ⟨$i|∂/∂R|$j⟩ - !!!!")
                            else
                                warn("Skipped smoothing of the slow oscilating curve ⟨$i|∂/∂R|$j⟩ found R₀ = $R₀ and τ₀ = $τ₀, ∫⟨$i|∂/∂R|$j⟩dR = $It, [$R¹, $Rᵉ]")
                                warn("Spectrum: $(abs(ddr_spectrum))")
                            end
                        end
                    end
                end
            end
            info("******************************************************")
        end
    end
    info("==========================================================================================")

    # ddr_spectrum = fft(ddr_sample)
    # Lˢ = length(ddr_spectrum)
    # ix_cutting = round(Int, Lˢ/80)
    # ix_cutting = ix_cutting ≥ 1 ? ix_cutting : 1
    # if ix_cutting <= 2
    #     warn("Possible rough smoothing with a single harmonic for ⟨$i|∂/∂R|$j⟩: [$ix_cutting, $(length(ddr_spectrum))]")
    # else
    #     info("Harmonics cutted: [$ix_cutting, $(length(ddr_spectrum))]")
    # end
    # ddr_spectrum[ix_cutting:end] = 0.0
    # ddr_sample = real(ifft(ddr_spectrum))
    # ∂_∂Rᴰᵈᵃᵗᵃ[l¹:lᵉ, k] = ddr_sample

    smoothed_∂_∂Rᵈ = matlupperx_ddr2matl(∂_∂Rᴰᵈᵃᵗᵃ)
    @assert size(∂_∂Rᵈ, 1) == size(smoothed_∂_∂Rᵈ, 1)
    @assert size(∂_∂Rᵈ, 2) == size(smoothed_∂_∂Rᵈ, 2)
    info("Smoothed matrix ∂_∂Rᵈ length = $(length(smoothed_∂_∂Rᵈ))")
    info("==========")

    info("Done.")
    return Rᶜ, Sᵛᵉᶜ, Hᵈ, smoothed_∂_∂Rᵈ, ∂_∂Rᵐ
end

function diabatize(Hₐ::Array{Function, 2}, ∂_∂R::Array{Function, 2}, ∂_∂Rᵐᵒᵈᵉˡ::Array{Function, 2},
  Rᵖᵒⁱⁿᵗˢ::Vector{Float64}, invert_R::Bool, Sˡ::Vector{Array{Float64, 2}}, use_prev_S_from::Nullable{Float64})
  Logging.configure(level=INFO)

  #increasing_order = Rᵖᵒⁱⁿᵗˢ[1] < Rᵖᵒⁱⁿᵗˢ[end]

  Nᵖᵒⁱⁿᵗˢ = size(Rᵖᵒⁱⁿᵗˢ, 1)
  Hᵈ = Vector{Array{Float64, 2}}(Nᵖᵒⁱⁿᵗˢ)
  ∂_∂Rᵈ = Vector{Array{Float64, 2}}(Nᵖᵒⁱⁿᵗˢ)
  Sᵛᵉᶜ = Vector{Array{Float64, 2}}(Nᵖᵒⁱⁿᵗˢ)

  Sᶠᵘⁿᶜ, S_spline = matl2matfsl(Rᵖᵒⁱⁿᵗˢ, Sˡ)
  Rᵛᵉᶜ = invert_R ? Rᵖᵒⁱⁿᵗˢ[end:-1:1] : Rᵖᵒⁱⁿᵗˢ
  Sˡᵛᵉᶜ = invert_R ? Sˡ[end:-1:1] : Sˡ
  N = size(Sˡᵛᵉᶜ[1], 1)
  Sᵖʳᵉᵛ = Array{Float64, 2}(N, N)
  info("Transforming matrix elements <|Ĥ|> and <|∂/∂R|> in interval [$(Rᵛᵉᶜ[1]), $(Rᵛᵉᶜ[end])]")
  for i = 1:Nᵖᵒⁱⁿᵗˢ
    R = Rᵛᵉᶜ[i]
    # ----
    use_prev_solution = !(isnull(use_prev_S_from) || (R > get(use_prev_S_from)))
    S = isnull(use_prev_S_from) || (R > get(use_prev_S_from)) ? Sˡᵛᵉᶜ[i] : Sᵖʳᵉᵛ
    if use_prev_solution
      info("Using previous transformation matrix at $R")
      S = round(S, 0)
    end
    S⁻¹ = S'
    # ----
    ∇S = Dierckx.derivative.(S_spline, R; nu=1)
    #∇S = matDerivative(R, S_spline)
    #∇S = Calculus.derivative.(Sᶠᵘⁿᶜ, R) # Calculus, what the f***???!!
    #∇S = dirtyDerivative.(Sᶠᵘⁿᶜ, R, 1e-6)
    Hᴬ = matf2mat(R, Hₐ); ∂_∂Rᴬ = matf2mat(R, ∂_∂R); ∂_∂Rᴹ = matf2mat(R, ∂_∂Rᵐᵒᵈᵉˡ)

    S⁻¹ = S'
    #S⁻¹ = inv(S)
    Hᴰ = S⁻¹*Hᴬ*S
    #∂_∂Rᴰ = S⁻¹*∂_∂Rᴬ*S + S⁻¹*∇S
    ∂_∂Rᴰ =
        if use_prev_solution
            ∇S .= 0.0 # We assume we have constants after some R,
            S⁻¹*∂_∂Rᴬ*S + S⁻¹* ∇S
        else
            S⁻¹ * (∂_∂Rᴬ - ∂_∂Rᴹ) * S
        end# + inv(S) * ∇S
    #∂_∂Rᴰ = ∂_∂Rᴬ

    Sᵛᵉᶜ[i] = S
    Hᵈ[i] = Hᴰ; ∂_∂Rᵈ[i] = ∂_∂Rᴰ

    Sᵖʳᵉᵛ = isnull(use_prev_S_from) || R > get(use_prev_S_from) ? S : Sᵖʳᵉᵛ
  end
  if invert_R
    return Rᵛᵉᶜ[end:-1:1], Hᵈ[end:-1:1], ∂_∂Rᵈ[end:-1:1], Sᵛᵉᶜ[end:-1:1]
  else
    return Rᵛᵉᶜ, Hᵈ, ∂_∂Rᵈ, Sᵛᵉᶜ
  end
end

function diabatize(Hₐ::Array{Function, 2}, ∂_∂R::Array{Function, 2},
  Rᵖᵒⁱⁿᵗˢ::Vector{Float64}, invert_R::Bool, Sˡ::Vector{Array{Float64, 2}}, use_prev_S_from::Nullable{Float64})
  Logging.configure(level=INFO)
  info("Diabatization with a precomputed transformation matrix")

  #increasing_order = Rᵖᵒⁱⁿᵗˢ[1] < Rᵖᵒⁱⁿᵗˢ[end]

  Nᵖᵒⁱⁿᵗˢ = size(Rᵖᵒⁱⁿᵗˢ, 1)
  Hᵈ = Vector{Array{Float64, 2}}(Nᵖᵒⁱⁿᵗˢ)
  ∂_∂Rᵈ = Vector{Array{Float64, 2}}(Nᵖᵒⁱⁿᵗˢ)
  Sᵛᵉᶜ = Vector{Array{Float64, 2}}(Nᵖᵒⁱⁿᵗˢ)

  Sᶠᵘⁿᶜ, S_spline = matl2matfsl(Rᵖᵒⁱⁿᵗˢ, Sˡ)
  Rᵛᵉᶜ = invert_R ? Rᵖᵒⁱⁿᵗˢ[end:-1:1] : Rᵖᵒⁱⁿᵗˢ
  Sˡᵛᵉᶜ = invert_R ? Sˡ[end:-1:1] : Sˡ
  N = size(Sˡᵛᵉᶜ[1], 1)
  Sᵖʳᵉᵛ = Array{Float64, 2}(N, N)
  info("Transforming matrix elements <|Ĥ|> and <|∂/∂R|> in interval [$(Rᵛᵉᶜ[1]), $(Rᵛᵉᶜ[end])]")
  progress = Progress(Nᵖᵒⁱⁿᵗˢ)
  for i = 1:Nᵖᵒⁱⁿᵗˢ
    R = Rᵛᵉᶜ[i]
    # ----
    S = isnull(use_prev_S_from) || (R > get(use_prev_S_from)) ? Sˡᵛᵉᶜ[i] : Sᵖʳᵉᵛ
    S⁻¹ = S'
    if !(isnull(use_prev_S_from) || (R > get(use_prev_S_from)))
      info("Using previous transformation matrix at $R")
    end
    # ----
    ∇S = Dierckx.derivative.(S_spline, R; nu=1)
    #∇S = matDerivative(R, S_spline)
    #∇S = Calculus.derivative.(Sᶠᵘⁿᶜ, R) # Calculus, what the f***???!!
    #∇S = dirtyDerivative.(Sᶠᵘⁿᶜ, R, 1e-6)
    Hᴬ = matf2mat(R, Hₐ); ∂_∂Rᴬ = matf2mat(R, ∂_∂R)

    Hᴰ = S⁻¹*Hᴬ*S
    ∂_∂Rᴰ = S⁻¹*∂_∂Rᴬ*S + S⁻¹*∇S
    #∂_∂Rᴰ = ∂_∂Rᴬ

    Sᵛᵉᶜ[i] = S
    Hᵈ[i] = Hᴰ; ∂_∂Rᵈ[i] = ∂_∂Rᴰ

    Sᵖʳᵉᵛ = isnull(use_prev_S_from) || R > get(use_prev_S_from) ? S : Sᵖʳᵉᵛ
    #info("Diabatization performed at the distance R = $R Bohr")
    ProgressMeter.next!(progress; showvalues = [(:index, i), (:distance, "$R, Bohr")])
  end
  if invert_R
    return Rᵛᵉᶜ[end:-1:1], Hᵈ[end:-1:1], ∂_∂Rᵈ[end:-1:1], Sᵛᵉᶜ[end:-1:1]
  else
    return Rᵛᵉᶜ, Hᵈ, ∂_∂Rᵈ, Sᵛᵉᶜ
  end
end

function transformationMatrix(Hₐ::Array{Function, 2},
  ∂_∂R::Array{Function, 2}, ∂_∂Rᵐᵒᵈᵉˡ::Array{Function, 2},
  Rᵛ::Vector{Float64},
  S₀ᵒʷⁿ::Nullable{Array{Float64, 2}},
  C::DiabatizationSettings)

  Logging.configure(level=INFO)

  # -----------
  N = size(Hₐ, 1)
  ℲR⃛ᵒᵖᵗ=(R₀::Float64, R::Float64, m₀, M₀, ΔRₘᵢₙ, ΔRₘₐₓ)->begin
    n = 1
    for i = 1:N, j = 1:N
      if i < j && j - i == 1
        n = max(splitn(R₀, R, ∂_∂R[i, j], m₀, M₀, ΔRₘᵢₙ, ΔRₘₐₓ), n)
      end
    end
    return splitxn(R₀, R, n)
  end
  # -----------

  Rˢᵗᵃʳᵗ = C.coordinate_start
  Rˢᵗᵒᵖ = C.coordinate_stop
  σR = sign(Rˢᵗᵒᵖ - Rˢᵗᵃʳᵗ)
  ΔRₘᵢₙ = σR * min(abs(C.coordinate_step[1]), abs(C.coordinate_step[2]))
  ΔRₘₐₓ = σR * max(abs(C.coordinate_step[1]), abs(C.coordinate_step[2]))

  # -----------
  progress = progressCreate(100, "Making R grid: ", :yellow); progressᶜ = 0
  Rᵖᵒⁱⁿᵗˢ = Vector{Float64}()
  Rₘᵢₙ = min(Rˢᵗᵃʳᵗ, Rˢᵗᵒᵖ); Rₘₐₓ = max(Rˢᵗᵃʳᵗ, Rˢᵗᵒᵖ)
  R⃜ = filter(R -> Rₘᵢₙ<=R<=Rₘₐₓ, σR > 0 ? Rᵛ : reverse(unique(Rᵛ))); L = length(R⃜)
  for l = 1:L-1
    Rˡ = R⃜[l]; Rˡ⁺¹ = R⃜[l + 1]
    R⃛ˡ = ℲR⃛ᵒᵖᵗ(Rˡ, Rˡ⁺¹, 2e-2, 1.2e4, ΔRₘᵢₙ, ΔRₘₐₓ)
    for k = 1:length(R⃛ˡ) - 1
      Rˡₖ = R⃛ˡ[k]
      push!(Rᵖᵒⁱⁿᵗˢ, Rˡₖ)
      progressᶜ = floor(Int, (1 - floor(abs(Rˡₖ - Rˢᵗᵒᵖ)/Rₘₐₓ)) * 100)
      progress!(progress, progressᶜ, [(:R, Rˡₖ), (:ΔR, abs(Rˡₖ - R⃛ˡ[k + 1]))])
    end
  end
  push!(Rᵖᵒⁱⁿᵗˢ, R⃜[L])
  finish!(progress)
  # -----------



  S₀ = isnull(S₀ᵒʷⁿ) ? eye(N, N) : get(S₀ᵒʷⁿ)
  info("Going to solve a Cauchy probem with initial conditions:\n$(S₀);\ncustom conditions are $(isnull(S₀ᵒʷⁿ) ? "null" : "not null")")
  S, Sᵈᵃᵗᵃ = problemCauchy(
    Rᵖᵒⁱⁿᵗˢ, S₀;
    prod_function = diabatizationODE_function,
    data = (N, ∂_∂R, ∂_∂Rᵐᵒᵈᵉˡ),
    ϵʳᵉˡ = 1e-5, ϵᵃᵇˢ = 1e-10
  )

  return Rᵖᵒⁱⁿᵗˢ[end:-1:1], S[end:-1:1], Sᵈᵃᵗᵃ[end:-1:1, 1:1:end]
end
