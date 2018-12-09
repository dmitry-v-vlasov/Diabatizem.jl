using Calculus
using Logging
using ProgressMeter
using Nullables

import Dierckx
import Sundials

abstract type AreaSystem end

mutable struct AreaBunch <: AreaSystem
    states::Vector{Int}
    areas::Vector{SinglePeakNonadiabaticArea}
    Rˡ::Float64
    Rʳ::Float64
end
function show(io::IO, ab::AreaBunch)
    print(io, "𝔹{⟨$(ab.states)⟩; Nᴬ = ($(length(ab.areas))); [$(ab.Rˡ), $(ab.Rʳ)]; $(ab.areas)}")
end

function solverTransformationMatrixForAreas(
    areas::Array{Vector{SinglePeakNonadiabaticArea}, 2},
    ∂_∂R::Array{Function, 2}, ∂_∂Rᵐᵒᵈᵉˡ::Array{Function, 2},
    Rᵛ::Vector{Float64}, C::DiabatizationSettings)
    Logging.configure(level=INFO)
    info("==== Starting selective transformation matrix solving ===")
    given_areas =
        sort(
            collect(
                Base.flatten(
                    filter(
                        vAⁱʲ -> !isempty(vAⁱʲ) && vAⁱʲ[1].states[1] < vAⁱʲ[1].states[2], areas)));
            by = A -> A.coordinate_∂_∂R, rev = true)
    info("The ordered area list:\n------")
    for area ∈ given_areas
        info("Area: $area")
    end
    info("\n------")
    sel_area_extras = Dict{Tuple{Int, Int, Float64}, Tuple{Float64, Float64}}()
    sel_area_bunch_exclude = Dict{Tuple{Int, Int, Float64}, Bool}()
    for sel_area ∈ C.areas
        s¹ = sel_area.states[1]; s² = sel_area.states[2]
        sel_area_key = (s¹, s², floor(sel_area.coordinate, 1))
        sel_area_extras[sel_area_key] = sel_area.extra_length
        sel_area_bunch_exclude[sel_area_key] = sel_area.bunch_exclude
    end
    info("Selected area extra lenghts:\n------")
    info("$sel_area_extras")
    info("\n------")
    ϵᵟᴿ = C.area_closeness
    info("Grouping the selected areas into bunches by their closeness within ϵᵟᴿ = $(ϵᵟᴿ) a.u.l.")
    bunched_areas = Set{Int}()
    bunches = Vector{AreaBunch}()
    for k = 1:length(given_areas)
        if k ∈ bunched_areas
            continue
        end
        Aᵏ = given_areas[k]
        Rᵏ⁻ˡ = Aᵏ.coordinate_from; Rᵏ⁻ʳ = Aᵏ.coordinate_to
        @assert(Rᵏ⁻ˡ < Rᵏ⁻ʳ, "Violated Rᵏ⁻ˡ < Rᵏ⁻ʳ: $Rᵏ⁻ˡ < $Rᵏ⁻ʳ")
        states_Aᵏ = Aᵏ.states
        extras_key_k = (states_Aᵏ[1], states_Aᵏ[2], floor(Aᵏ.coordinate_∂_∂R, 1))
        extras_Aᵏ = sel_area_extras[extras_key_k]
        Rᵏ⁻ˡᵉ = Rᵏ⁻ˡ - extras_Aᵏ[1]; Rᵏ⁻ʳᵉ = Rᵏ⁻ʳ + extras_Aᵏ[2];
        #info("R-Extras for $extras_key_k ⇒ $extras_Aᵏ; Interval becomes [$Rᵏ⁻ˡᵉ, $Rᵏ⁻ʳᵉ] ")

        remaining_areas = view(given_areas, (k + 1):length(given_areas))
        remaining_bunch_of_areas =
            map(j_Aʲ -> j_Aʲ[2],
                filter(
                    j_Aʲ -> begin
                        info("Last bunched area: $extras_key_k, interval - [$Rᵏ⁻ˡᵉ |, $Rᵏ⁻ˡ, $Rᵏ⁻ʳ|, $Rᵏ⁻ʳᵉ]")

                        j, Aʲ = j_Aʲ
                        states_Aʲ = Aʲ.states
                        Rʲ⁻ˡ = Aʲ.coordinate_from; Rʲ⁻ʳ = Aʲ.coordinate_to
                        if states_Aʲ == states_Aᵏ || states_Aʲ[1] == states_Aᵏ[2] || states_Aʲ[2] == states_Aᵏ[1]
                            extras_key_j = (states_Aʲ[1], states_Aʲ[2], floor(Aʲ.coordinate_∂_∂R, 1))
                            extras_Aʲ = sel_area_extras[extras_key_j]
                            Rʲ⁻ˡᵉ = Rʲ⁻ˡ - extras_Aʲ[1]; Rʲ⁻ʳᵉ = Rʲ⁻ʳ + extras_Aʲ[2];

                            @assert(Rʲ⁻ˡ < Rʲ⁻ʳ, "Violated Rʲ⁻ˡ < Rʲ⁻ʳ: $Rʲ⁻ˡ < $Rʲ⁻ʳ")
                            b_j_k_overlapped = (Rᵏ⁻ʳ > Rʲ⁻ˡ && Rᵏ⁻ˡ < Rʲ⁻ʳ) || (Rʲ⁻ʳ > Rᵏ⁻ˡ && Rʲ⁻ˡ < Rᵏ⁻ʳ)
                            b_j_k_closed = (Rᵏ⁻ʳ > Rʲ⁻ˡ && Rᵏ⁻ˡ > Rʲ⁻ʳ && (Rᵏ⁻ˡ - Rʲ⁻ʳ) < ϵᵟᴿ) ||
                                        (Rʲ⁻ʳ > Rᵏ⁻ˡ && Rʲ⁻ˡ > Rᵏ⁻ʳ && (Rʲ⁻ʳ - Rᵏ⁻ˡ) < ϵᵟᴿ)
                            b_j_k_overlapped_e = (Rᵏ⁻ʳᵉ > Rʲ⁻ˡᵉ && Rᵏ⁻ˡᵉ < Rʲ⁻ʳᵉ) || (Rʲ⁻ʳᵉ > Rᵏ⁻ˡᵉ && Rʲ⁻ˡᵉ < Rᵏ⁻ʳᵉ)
                            b_j_k_closed_e = (Rᵏ⁻ʳᵉ > Rʲ⁻ˡᵉ && Rᵏ⁻ˡᵉ > Rʲ⁻ʳᵉ && (Rᵏ⁻ˡᵉ - Rʲ⁻ʳᵉ) < ϵᵟᴿ) ||
                                        (Rʲ⁻ʳᵉ > Rᵏ⁻ˡᵉ && Rʲ⁻ˡᵉ > Rᵏ⁻ʳᵉ && (Rʲ⁻ʳᵉ - Rᵏ⁻ˡᵉ) < ϵᵟᴿ)
                            #info("Test $states_Aᵏ <-> $states_Aʲ: o=$b_j_k_overlapped, c=$b_j_k_closed, $(Rᵏ⁻ˡ - Rʲ⁻ʳ), $ϵᵟᴿ, $Rᵏ⁻ˡ, $Rʲ⁻ʳ)")
                            if b_j_k_overlapped || b_j_k_closed
                                #info("Closed areas for states $states_Aᵏ: overlapped=$b_j_k_overlapped, closed=$b_j_k_closed; j:[$Rʲ⁻ˡ, $Rʲ⁻ʳ], k:[$Rᵏ⁻ˡ, $Rᵏ⁻ʳ]")
                                if sel_area_bunch_exclude[extras_key_j]
                                    warn("Bunch exclude with key $extras_key_j")
                                    return false
                                end
                                push!(bunched_areas, k + j)
                                info("The areas $extras_key_k and $extras_key_j are bunched (intervals)")
                            elseif b_j_k_overlapped_e || b_j_k_closed_e
                                #info("Extra-Closed areas for states $states_Aᵏ: overlapped=$b_j_k_overlapped, closed=$b_j_k_closed; j:[$Rʲ⁻ˡ, $Rʲ⁻ʳ], k:[$Rᵏ⁻ˡ, $Rᵏ⁻ʳ]")
                                #info("R-Extras for $extras_key_j ⇒ $extras_Aʲ; Interval becomes [$Rʲ⁻ˡᵉ, $Rʲ⁻ʳᵉ] ")
                                if sel_area_bunch_exclude[extras_key_j]
                                    warn("Bunch exclude with key $extras_key_j")
                                    return false
                                end
                                push!(bunched_areas, k + j)
                                info("The areas $extras_key_k and $extras_key_j are bunched (extra intervals)")
                            else
                                warn(
                                """
                                    The areas $extras_key_k and $extras_key_j are not bunched:
                                    $states_Aᵏ interval: [$Rᵏ⁻ˡᵉ |, $Rᵏ⁻ˡ, $Rᵏ⁻ʳ|, $Rᵏ⁻ʳᵉ]
                                    $states_Aʲ interval: [$Rʲ⁻ˡᵉ |, $Rʲ⁻ˡ, $Rʲ⁻ʳ|, $Rʲ⁻ʳᵉ]
                                """)
                            end
                            if b_j_k_overlapped || b_j_k_closed || b_j_k_overlapped_e || b_j_k_closed_e
                                info("Changing last bunched area $extras_key_k→$extras_key_j")
                                Rᵏ⁻ˡ = Aʲ.coordinate_from; Rᵏ⁻ʳ = Aʲ.coordinate_to
                                states_Aᵏ = states_Aʲ
                                extras_key_k = extras_key_j
                                extras_Aᵏ = sel_area_extras[extras_key_k]
                                Rᵏ⁻ˡᵉ = Rᵏ⁻ˡ - extras_Aᵏ[1]; Rᵏ⁻ʳᵉ = Rᵏ⁻ʳ + extras_Aᵏ[2];
                            end
                            return Bool(b_j_k_overlapped || b_j_k_closed || b_j_k_overlapped_e || b_j_k_closed_e)
                        else
                            return false
                        end
                    end,
                    enumerate(remaining_areas)))
        bunch_of_areas = Vector{SinglePeakNonadiabaticArea}()
        push!(bunch_of_areas, Aᵏ)
        for remaining_area_from_bunch ∈ remaining_bunch_of_areas
            push!(bunch_of_areas, remaining_area_from_bunch)
        end
        Rᵇ⁻ˡ = minimum(map(Aᵇ -> Aᵇ.coordinate_from, bunch_of_areas))
        Rᵇ⁻ʳ = maximum(map(Aᵇ -> Aᵇ.coordinate_to, bunch_of_areas))
        states_Aᵇ = sort(unique(Base.flatten(map(Aᵇ -> Aᵇ.states, bunch_of_areas))))
        new_bunch = AreaBunch(states_Aᵇ, bunch_of_areas, Rᵇ⁻ˡ, Rᵇ⁻ʳ)
        #info("New area bunch: $new_bunch")
        push!(bunches, new_bunch)
    end
    info("The bunched areas are the following:\n------")
    for bunch ∈ bunches
        info("B-Areas: $bunch")
    end
    info("\n------")
    info("Solving Cauchy problems...")

    last_matrices = Dict{Vector{Int}, Array{Float64, 2}}()

    solutions = Vector{LocalSolution}();
    for bunch ∈ bunches
        info("=====================")
        info("Making a new Cauchy problem for the areas:\n$bunch")
        states = bunch.states
        info("The involved states: ⟨$states⟩")
        extra_lengths = Vector{Tuple{Float64, Float64}}()
        foreach(
            Aᵇ -> begin
                extras_key_b = (Aᵇ.states[1], Aᵇ.states[2], floor(Aᵇ.coordinate_∂_∂R, 1))
                info("Trying key: $extras_key_b")
                extras_Aᵇ = sel_area_extras[extras_key_b]
                info("Key result: $extras_Aᵇ")
                push!(extra_lengths, extras_Aᵇ)
            end,
            sort(bunch.areas; by = Aᵇ -> Aᵇ.coordinate_∂_∂R))
        info("Result extra lenghts: $extra_lengths")
        extra_left = extra_lengths[1][1]
        extra_right = extra_lengths[end][2]
        info("The areas extra lengths: left = $extra_left, right = $extra_right")
        info("From-s: $(map(Aⁱ -> Aⁱ.coordinate_from, bunch.areas))")
        info("To-s: $(map(Aⁱ -> Aⁱ.coordinate_to, bunch.areas))")
        Rˢᵗᵃʳᵗ = bunch.Rʳ + extra_right      # From higher R
        Rˢᵗᵒᵖ = bunch.Rˡ - extra_left        # To lower R
        info("The solution interval: [$(Rˢᵗᵃʳᵗ), $(Rˢᵗᵒᵖ)] derived from [$(bunch.Rʳ), $(bunch.Rˡ)]")
        s¹ = states[1]; s² = states[end]
        Nˡᵒᶜ = s² - s¹ + 1
        view_∂_∂Rˡᵒᶜ = view(∂_∂R, s¹:s², s¹:s²)
        view_∂_∂Rᵐᵒᵈᵉˡ⁻ˡᵒᶜ = view(∂_∂Rᵐᵒᵈᵉˡ, s¹:s², s¹:s²)
        ∂_∂Rˡᵒᶜ = Array{Function, 2}(view_∂_∂Rˡᵒᶜ)
        ∂_∂Rᵐᵒᵈᵉˡ⁻ˡᵒᶜ = Array{Function, 2}(view_∂_∂Rᵐᵒᵈᵉˡ⁻ˡᵒᶜ)
        info("Check: $(∂_∂Rˡᵒᶜ)")
        info("The local solution size is $(Nˡᵒᶜ)×$(Nˡᵒᶜ).")
        S₀ˡᵒᶜ = eye(Nˡᵒᶜ, Nˡᵒᶜ)

        if C.keep_initial_conditions
            if haskey(last_matrices, states)
                S₀ˡᵒᶜ = last_matrices[states]
                info("For |$states⟩ - Found initial condition matrix from the previous area:\n$S₀ˡᵒᶜ")
            else
                #info("For |$states⟩ - Identity matrix is used:\n$S₀ˡᵒᶜ")
                if length(states) == 2
                    info("Detected a couple of neigbour solutions for the states |$states⟩.")
                    entries = collect(last_matrices)
                    ix_lower = findfirst(entry -> entry[1][end] == states[1], entries)
                    ix_upper = findfirst(entry -> entry[1][1] == states[2], entries)
                    if ix_upper ≠ 0
                        Sᵁᵖ = entries[ix_upper][2]
                        info("There exists a previous solution for the state |$(states[2])⟩: $(Sᵁᵖ[1, 1]) (upper)")
                        S₀ˡᵒᶜ[end, end] = Sᵁᵖ[1, 1]
                    end
                    if ix_lower ≠ 0
                        Sᴸᵒ = entries[ix_lower][2]
                        info("There exists a previous solution for the state |$(states[1])⟩: $(Sᴸᵒ[end, end]) (lower)")
                        S₀ˡᵒᶜ[1, 1] = Sᴸᵒ[end, end]
                    end
                    info("For |$states⟩ - The following matrix is used:\n$S₀ˡᵒᶜ")
                else
                    info("For |$states⟩ - Identity matrix is used:\n$S₀ˡᵒᶜ")
                end
            end
        else
        end

        info("Solving... (Nˡᵒᶜ = $(Nˡᵒᶜ))")
        Rᵖᵒⁱⁿᵗˢ, S, Sᵈᵃᵗᵃ = transformationMatrixForArea(Rˢᵗᵃʳᵗ, Rˢᵗᵒᵖ, ∂_∂Rˡᵒᶜ, ∂_∂Rᵐᵒᵈᵉˡ⁻ˡᵒᶜ, Rᵛ, Nullable(S₀ˡᵒᶜ), C)
        info("The solution contains $(length(Rᵖᵒⁱⁿᵗˢ)) points in the interval [$(Rᵖᵒⁱⁿᵗˢ[1]), $(Rᵖᵒⁱⁿᵗˢ[end])].")
        N = size(∂_∂Rᵐᵒᵈᵉˡ, 1)
        info("Extending the solution matrix size from $(Nˡᵒᶜ)×$(Nˡᵒᶜ) to $N×$N")
        ext_S = Vector{Array{Float64, 2}}()
        for l ∈ 1:length(Rᵖᵒⁱⁿᵗˢ)
            ext_Sˡ = eye(N, N)
            for i = s¹:s², j = s¹:s²
                ext_Sˡ[i, j] = S[l][i - s¹ + 1, j - s¹ + 1]
            end
            push!(ext_S, ext_Sˡ)
        end
        info("Done")
        R¹ = Rᵖᵒⁱⁿᵗˢ[1]; Rᵉ = Rᵖᵒⁱⁿᵗˢ[end]

        peaks = Vector{Tuple{Float64, Float64}}()
        for area ∈ bunch.areas
            R₀ᵢ = area.coordinate_∂_∂R
            τᴿ⁰ⁱ = area.value_∂_∂R
            peak = (R₀ᵢ, τᴿ⁰ⁱ)
            push!(peaks, peak)
        end
        info("Peaks in the local solution: $(peaks)")

        info("Adding to solution collection.")
        solution = LocalSolution(states,
            (Rᵖᵒⁱⁿᵗˢ[1], Rᵖᵒⁱⁿᵗˢ[end]),
            Rᵖᵒⁱⁿᵗˢ, Rᵖᵒⁱⁿᵗˢ,
            peaks,
            ext_S, Sᵈᵃᵗᵃ);
        last_matrix = Array{Float64, 2}(view(solution.S[1], s¹:s², s¹:s²))
        last_matrices[states] = round(last_matrix)
        info("For the states ⟨$states⟩ the last matrix at R = $(solution.points[1]) is\n$(mat2string(last_matrix))\nrounded to\n$(mat2string(last_matrices[states]))")
        push!(solutions, solution)
        info("Done")
    end
    info("End of making Cauchy problems")
    info("==== End of selective transformation matrix solving ===")
    return solutions;
end

function transformationMatrixForArea(
    Rˢᵗᵃʳᵗ::Float64, Rˢᵗᵒᵖ::Float64,
    ∂_∂R::Array{Function, 2}, ∂_∂Rᵐᵒᵈᵉˡ::Array{Function, 2},
    Rᵛ::Vector{Float64},
    S₀ᵒʷⁿ::Nullable{Array{Float64, 2}},
    C::DiabatizationSettings)

  Logging.configure(level=INFO)
  info("Transformation matrix calculation")

  # -----------
  N = size(∂_∂Rᵐᵒᵈᵉˡ, 1)
  ℲR⃛ᵒᵖᵗ = (R₀::Float64, R::Float64, m₀, M₀, ΔRₘᵢₙ, ΔRₘₐₓ) -> begin
    n = 1
    for i = 1:N, j = 1:N
      if i < j && j - i == 1
        nˡ = max(splitn(R₀, R, ∂_∂R[i, j], m₀, M₀, ΔRₘᵢₙ, ΔRₘₐₓ), n)
        if nˡ > n
            n = nˡ
        end
      end
    end
    return splitxn(R₀, R, n)
  end
  # -----------

  σR = sign(Rˢᵗᵒᵖ - Rˢᵗᵃʳᵗ)
  ΔRₘᵢₙ = σR * min(abs(C.coordinate_step[1]), abs(C.coordinate_step[2]))
  ΔRₘₐₓ = σR * max(abs(C.coordinate_step[1]), abs(C.coordinate_step[2]))

  # -----------
  progress = progressCreate(100, "Making R grid: ", :yellow); progressᶜ = 0
  Rᵖᵒⁱⁿᵗˢ = Vector{Float64}()
  Rₘᵢₙ = min(Rˢᵗᵃʳᵗ, Rˢᵗᵒᵖ); Rₘₐₓ = max(Rˢᵗᵃʳᵗ, Rˢᵗᵒᵖ)
  R⃜ = filter(R -> Rₘᵢₙ<=R<=Rₘₐₓ, σR > 0 ? unique(Rᵛ) : reverse(unique(Rᵛ))); L = length(R⃜)
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
  if(!allunique(Rᵖᵒⁱⁿᵗˢ))
      Rᵖᵒⁱⁿᵗˢ⁻ᵗ = unique(Rᵖᵒⁱⁿᵗˢ)Rᵖᵒⁱⁿᵗˢ[end:-1:1], S[end:-1:1], Sᵈᵃᵗᵃ[end:-1:1, 1:1:end]
      warn("Duplicate points are found. Resizing their vector from $(length(Rᵖᵒⁱⁿᵗˢ)) to $(length(Rᵖᵒⁱⁿᵗˢ⁻ᵗ)).")
      resize!(Rᵖᵒⁱⁿᵗˢ, length(Rᵖᵒⁱⁿᵗˢ⁻ᵗ))
      Rᵖᵒⁱⁿᵗˢ = Rᵖᵒⁱⁿᵗˢ⁻ᵗ
  end
  finish!(progress)
  # -----------


  S₀ = isnull(S₀ᵒʷⁿ) ? eye(N, N) : get(S₀ᵒʷⁿ)
  info("Going to solve a Cauchy problem with initial conditions:\n$(mat2string(S₀));\ncustom conditions are $(isnull(S₀ᵒʷⁿ) ? "null" : "not null")")
  info("The size N = $N")
  S, Sᵈᵃᵗᵃ = problemCauchy(
    Rᵖᵒⁱⁿᵗˢ, S₀;
    prod_function = diabatizationODE_function,
    data = (N, ∂_∂R, ∂_∂Rᵐᵒᵈᵉˡ),
    ϵʳᵉˡ = 1e-10, ϵᵃᵇˢ = 1e-12
  )

  if σR > 0
      info("Direct result")
      return Rᵖᵒⁱⁿᵗˢ, S, Sᵈᵃᵗᵃ
  else
      info("Inverted result")
      return Rᵖᵒⁱⁿᵗˢ[end:-1:1], S[end:-1:1], Sᵈᵃᵗᵃ[end:-1:1, 1:1:end]
  end
end

function error_S(S::Vector{Array{Float64, 2}})
  L = size(S, 1)
  N = size(S[1], 1)
  ϵ_S = Vector{Array{Float64, 2}}(L)
  for l = 1:L
    ϵ_S[l] = S[l]'*S[l]
  end
  ϵ_Sᵈᵃᵗᵃ = matl2mdata(ϵ_S)
  return ϵ_S, ϵ_Sᵈᵃᵗᵃ
end

let
  global diabatizationODE_function
  S = Array{Float64, 2}(); dS_dR = Array{Float64, 2}()
  ∂_∂R = Array{Float64, 2}(); ∂_∂Rᵐᵒᵈᵉˡ = Array{Float64, 2}()

  """
  We do not aspite to have an optimal implementation and
  we prefer to have this function with vivid matrix formulae.
  """
  function diabatizationODE_function(R::Float64, S_v::Vector{Float64}, dS_dR_v::Vector{Float64}, data::Tuple{Int, Array{Function, 2}, Array{Function, 2}})
    F_∂_∂R = data[2]; F_∂_∂Rᵐᵒᵈᵉˡ = data[3]

    N = data[1]
    @assert(N*N == size(S_v, 1), "With N = $N: $(N*N) ≠ $(size(S_v, 1)).")
    if size(S, 1) ≠ N || size(dS_dR, 1) ≠ N || size(∂_∂R, 1) ≠ N || size(∂_∂Rᵐᵒᵈᵉˡ, 1) ≠ N
      S = Array{Float64, 2}(N, N); dS_dR = Array{Float64, 2}(N, N)
      ∂_∂R = Array{Float64, 2}(N, N); ∂_∂Rᵐᵒᵈᵉˡ = Array{Float64, 2}(N, N)
    end

    vec2mat!(S_v, S)
    for i = 1:N, j = 1:N
      ∂_∂R[i, j] = F_∂_∂R[i, j](R)
      ∂_∂Rᵐᵒᵈᵉˡ[i, j] = F_∂_∂Rᵐᵒᵈᵉˡ[i, j](R)
    end

    #dS_dR = (S*∂_∂R - ∂_∂R*S) - S*∂_∂Rᵐᵒᵈᵉˡ
    dS_dR = -∂_∂Rᵐᵒᵈᵉˡ*S

    mat2vec!(dS_dR, dS_dR_v)

    return Sundials.CV_SUCCESS
  end
end

function problemCauchy(
  Xᵖᵒⁱⁿᵗˢ::Vector{Float64},
  Y₀::Array{Float64, 2};
  prod_function::Function = nothing,
  data::Tuple{Int, Array{Function, 2}, Array{Function, 2}} = nothing,
  ϵʳᵉˡ::Float64 = 1e-10,
  ϵᵃᵇˢ::Float64 = 1e-12)

  N = size(Y₀, 1)

  Yⁱⁿⁱᵗ = Vector{Float64}(N*N); fill!(Yⁱⁿⁱᵗ, 0)
  for i = 1:N, j = 1:N
    l = mvec(i, j, N); Yⁱⁿⁱᵗ[l] = Y₀[i, j]
  end

  Yʳᵉˢ = Sundials.cvode(
    (Xᵖ, Yᵛ, dYᵛ_dx) -> prod_function(Xᵖ, Yᵛ, dYᵛ_dx, data),
    Yⁱⁿⁱᵗ, Xᵖᵒⁱⁿᵗˢ, nothing;
    integrator = :BDF, reltol = ϵʳᵉˡ, abstol = ϵᵃᵇˢ)

  Nᵖᵒⁱⁿᵗˢ = size(Yʳᵉˢ, 1)
  Nˢᵒˡᵘᵗⁱᵒⁿˢ = size(Yʳᵉˢ, 2)
  if Nˢᵒˡᵘᵗⁱᵒⁿˢ ≠ N*N throw(ErrorException("Unexpected: Nˢᵒˡᵘᵗⁱᵒⁿˢ ≠ N×N: $Nˢᵒˡᵘᵗⁱᵒⁿˢ ≠ $N×$N ($(N*N))")) end
  Y = Vector{Array{Float64, 2}}(Nᵖᵒⁱⁿᵗˢ)
  for k = 1:Nᵖᵒⁱⁿᵗˢ
    Yₖ = Array{Float64, 2}(N, N)
    for l = 1:N*N
      i, j = mpos(l, N)
      Yₖ[i, j] = Yʳᵉˢ[k, l]
    end
    Y[k] = Yₖ
  end

  return Y, Yʳᵉˢ
end
