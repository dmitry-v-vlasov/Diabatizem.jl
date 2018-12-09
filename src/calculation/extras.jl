
import Dierckx

function expandLocalSolutions(solutions::Vector{LocalSolution}, Rᵛ::Vector{Float64}, C::DiabatizationSettings)
    for solution ∈ solutions
        @info "Working with solution:\n$solution"
        S¹ = solution.S[1]; Sᵉ = solution.S[end]
        R¹ = solution.points[1]; Rᵉ = solution.points[end]
        S¹ᵗ = round(S¹); Sᵉᵗ = round(Sᵉ)
        Rᵛ¹ = Rᵛ[1]; Rᵛᵉ = Rᵛ[end]
        Svf = matl2matfsl(solution.points, solution.S)[1]
        @info "=========================================================================================="
        @info "******************************************************************************************"
        @info "Boundary matrices:\n - S¹=\n$(mat2string(S¹));\n - Sᵉ=\n$(mat2string(Sᵉ))"
        @info "-----"
        @info "Argument boundaries: [$R¹, $Rᵉ]"
        @info "+++++"
        @info "Target boundary matrices:\n - S¹ᵗ=\n$(mat2string(S¹ᵗ));\n - Sᵉᵗ=\n$(mat2string(Sᵉᵗ))"
        @info "-----"
        @info "Ταrget argument boundaries: [$Rᵛ¹, $Rᵛᵉ]"
        @info "******************************************************************************************"
        Rᴸˢ, Rᴿˢ, Sᴸˢ, Sᴿˢ = smoothing_sigmoid(solution, Rᵛ, S¹ᵗ, Sᵉᵗ, 5e-2, 5e-2)
        @info "Expanding the solution from the interval [$(Rᴸˢ[1]), $(Rᴿˢ[end])] to [$(Rᵛ[1]), $(Rᵛ[end])]"
        Rᵐ = mergeGrids(Rᵛ, Rᴸˢ, solution.points, Rᴿˢ)
        @assert issorted(Rᵐ)
        Rᵖ = solution.points
        Rˢ = Vector{Float64}(undef, 0)
        Sˢ = Vector{Matrix{Float64}}(undef, 0)
        for R ∈ Rᵐ
            push!(Rˢ, R)
            if R < Rᴸˢ[1]
                push!(Sˢ, S¹ᵗ)
            elseif Rᴸˢ[1] ≤ R ≤ Rᴸˢ[end]
                push!(Sˢ, matf2mat(R, Sᴸˢ))
            elseif Rᴸˢ[end] < R < Rᴿˢ[1]
                push!(Sˢ, matf2mat(R, Svf))
            elseif Rᴿˢ[1] ≤ R ≤ Rᴿˢ[end]
                push!(Sˢ, matf2mat(R, Sᴿˢ))
            elseif R > Rᴿˢ[end]
                push!(Sˢ, Sᵉᵗ)
            end
        end
        solution.points = Rˢ
        solution.S = Sˢ
        @info "The result solution size is $(length(solution.S)) in with $(length(solution.points)) points"
        @info "=========================================================================================="
    end
end

function mergeGrids(Rᵛ, Rᴸˢ, Rˢ, Rᴿˢ)
    Rᵐ = Vector{Float64}(undef, 0)
    @assert Rᵛ[1] <= Rᴸˢ[1]
    @assert Rᴸˢ[1] < Rˢ[1]
    @assert Rˢ[1] < Rˢ[end]
    @assert Rˢ[end] < Rᴿˢ[end]
    @assert Rᴿˢ[end] <= Rᵛ[end]
    # 1. left interval
    iR = findlast(R -> R < Rᴸˢ[1], Rᵛ)
    if iR > 0
        @info "Merging left interval [$(Rᵛ[1]), $(Rᵛ[iR])]"
        for i = 1:iR
            push!(Rᵐ, Rᵛ[i])
        end
    end
    # 2. left boundary
    @info "Merging left boundary [$(Rᴸˢ[1]), $(Rᴸˢ[end])]"
    for Rᴸˢⁱ ∈ Rᴸˢ
        push!(Rᵐ, Rᴸˢⁱ)
    end
    # 3. main interval
    iRᴸ = findfirst(R -> R > Rᴸˢ[end], Rˢ)
    iRᴿ = findlast(R -> R < Rᴿˢ[1], Rˢ)
    @assert iRᴸ > 0
    @assert iRᴿ > 0
    @assert iRᴸ < iRᴿ
    @info "Merging solution interval [$(Rˢ[iRᴸ]), $(Rˢ[iRᴿ])]"
    for i = iRᴸ:iRᴿ
        push!(Rᵐ, Rˢ[i])
    end
    # 4. right boundary
    @info "Merging right boundary [$(Rᴿˢ[1]), $(Rᴿˢ[end])]"
    for Rᴿˢᵢ ∈ Rᴿˢ
        push!(Rᵐ, Rᴿˢᵢ)
    end
    # 5. right interval
    iR = findfirst(R -> R > Rᴿˢ[end], Rᵛ)
    if iR > 0
        @info "Merging right interval [$(Rᵛ[iR]), $(Rᵛ[end])]"
        for i = iR:length(Rᵛ)
            push!(Rᵐ, Rᵛ[i])
        end
    end
    # ----
    return unique(Rᵐ)
end

function smoothing_sigmoid(
    solution::LocalSolution, Rᵛ::Vector{Float64},
    S¹ᵗ::Matrix{Float64}, Sᵉᵗ::Matrix{Float64},
    ϵ::Float64, ϵ⁽¹⁾::Float64)
    states = solution.states; s¹ = states[1]; sᵉ = states[end]
    @info "ΣSigmoid smoothingΣ for $states"
    R = solution.points

    Sv = solution.S
    S = Vector{Matrix{Float64}}(undef, 0)
    for Si ∈ Sv
        push!(S, Matrix{Float64}(view(Si, s¹:sᵉ, s¹:sᵉ)))
    end
    Svf = matl2matfsl(R, Sv)[1]

    Rᴸⁱⁿᵗ = solutionBoundaryPoints(R, S, ϵ, ϵ⁽¹⁾)
    Rᴿⁱⁿᵗ = solutionBoundaryPoints(R, S, 0.15, 0.15; rev=true, max_Δ⁽⁰⁾ = 40, max_Δ⁽¹⁾ = 40)
    Rᴸᵇ = Rᴸⁱⁿᵗ[1]
    Rᴿᵇ = Rᴿⁱⁿᵗ[end]
    Rᴸ, Rᴿ = argumentGrid(Rᴸⁱⁿᵗ, Rᴿⁱⁿᵗ, Rᵛ)
    @assert issorted(Rᴸ)
    @assert issorted(Rᴿ)
    Rᴸᵉˣᵗ = filter(R -> R < Rᴸᵇ, Rᴸ)
    Rᴿᵉˣᵗ = filter(R -> R > Rᴿᵇ, Rᴿ)

    @info "The resulting boundary intervals: left - [$(Rᴸ[1]), $(Rᴸ[end])], right - [$(Rᴿ[1]), $(Rᴿ[end])]"
    Sᴸⁱⁿᵗ⁻ᶠ, Sᴿⁱⁿᵗ⁻ᶠ = internalFunctions(Rᴸⁱⁿᵗ, Rᴿⁱⁿᵗ, Sv)
    Sᴸᵛ, Sᴿᵛ = valuableFunctions(Rᴸⁱⁿᵗ, Rᴿⁱⁿᵗ, Sᴸⁱⁿᵗ⁻ᶠ, Sᴿⁱⁿᵗ⁻ᶠ, Rᴸ, Rᴿ, s¹, sᵉ)
    Sᴸᵈ, Sᴿᵈ = dummyFunctions(S¹ᵗ, Sᵉᵗ, size(S¹ᵗ, 1))

    @info "Making boundary sigmoid functions."
    @info "Left interval: [$(Rᴸ[1]) ..[$(length(Rᴸᵉˣᵗ))].., |$(Rᴸᵇ) ...[$(length(Rᴸⁱⁿᵗ))]... $(Rᴸ[end])]"
    @info "Right interval: [$(Rᴿ[1]) ...[$(length(Rᴿⁱⁿᵗ))]... $(Rᴿᵇ)|, ..[$(length(Rᴿᵉˣᵗ))].. $(Rᴿ[end])]"
    Sᴸ, Sᴿ = boundaryFunctions(Rᴸ, Rᴿ, Sᴸⁱⁿᵗ⁻ᶠ, Sᴿⁱⁿᵗ⁻ᶠ, Sᴸᵈ, Sᴿᵈ)
    @info "The boundary functions are ready."
    @info "++++++"
    @info "Sᴸ(R=$(Rᴸ[1])) =\n$(mat2string(matf2mat(Rᴸ[1], Sᴸ)))"
    @info "---"
    @info "Sᴿ(R=$(Rᴿ[end])) =\n$(mat2string(matf2mat(Rᴿ[end], Sᴿ)))"
    @info "++++++"

    @info "----------"
    # αˢʰᵃʳᵖ = 4.0
    # x₀ = (X[1] + X[end])/2; α = σf * abs(X[end]-X[1]) / αˢʰᵃʳᵖ
    @info "----------"
    return Rᴸ, Rᴿ, Sᴸ, Sᴿ
end

function boundaryFunctions(
    Rᴸ::Vector{Float64}, Rᴿ::Vector{Float64},
    Sᴸᵛ::Array{Function, 2}, Sᴿᵛ::Array{Function, 2},
    Sᴸᵈ::Array{Function, 2}, Sᴿᵈ::Array{Function, 2})
    @info "Making boundary functions."
    N = size(Sᴸᵛ, 1);
    Sᴸ = Array{Function, 2}(N, N); Sᴿ = Array{Function, 2}(N, N)

    for i = 1:N, j = 1:N
        # ----
        αˢʰᵃʳᵖ = 10.0
        σf = sign(Sᴸᵈ[i, j](Rᴸ[1]) - Sᴸᵛ[i, j](Rᴸ[end]))
        R₀ = inv_ϕ₀ * Rᴸ[1] + ϕ₀ * Rᴸ[end]
        α = abs(Rᴸ[end] - Rᴸ[1]) / αˢʰᵃʳᵖ
        # ----
        Sᴸ[i, j] = sigmoid_of_name(Sᴸᵈ[i, j], Sᴸᵛ[i, j], R₀, α, "Sᴸ[$i, $j]")
    end

    for i = 1:N, j = 1:N
        # ----
        αˢʰᵃʳᵖ = 10.0
        σf = sign(Sᴿᵛ[i, j](Rᴿ[end]) - Sᴿᵈ[i, j](Rᴿ[1]))
        R₀ = ϕ₀ * Rᴿ[1] + inv_ϕ₀ * Rᴿ[end]
        α = abs(Rᴿ[end] - Rᴿ[1]) / αˢʰᵃʳᵖ
        # ----
        Sᴿ[i, j] = sigmoid_of_name(Sᴿᵛ[i, j], Sᴿᵈ[i, j], R₀, α, "Sᴿ[$i, $j]")
    end

    return Sᴸ, Sᴿ
end

function dummyFunctions(
    S¹ᵗ::Matrix{Float64}, Sᵉᵗ::Matrix{Float64}, N::Int)
    Sᴸ = Array{Function, 2}(N, N); Sᴿ = Array{Function, 2}(N, N)
    for i = 1:N, j = 1:N
        Sᴸ[i, j] = R -> S¹ᵗ[i, j]
        Sᴿ[i, j] = R -> Sᵉᵗ[i, j]
    end
    return Sᴸ, Sᴿ
end

function valuableFunctions(
    Rᴸⁱⁿᵗ::Vector{Float64}, Rᴿⁱⁿᵗ::Vector{Float64},
    Sᴸⁱⁿᵗ⁻ᶠ::Array{Function, 2}, Sᴿⁱⁿᵗ⁻ᶠ::Array{Function, 2},
    Rᴸ::Vector{Float64}, Rᴿ::Vector{Float64},
    s¹::Int, sᵉ::Int)
    @info "Making valuable functions."
    @assert size(Sᴸⁱⁿᵗ⁻ᶠ, 1) == size(Sᴿⁱⁿᵗ⁻ᶠ, 1)
    N = size(Sᴸⁱⁿᵗ⁻ᶠ, 1);
    Rᴸ⁰ = Rᴸⁱⁿᵗ[1]; Rᴸᵉ = Rᴸⁱⁿᵗ[end]
    Rᴿ⁰ = Rᴿⁱⁿᵗ[1]; Rᴿᵉ = Rᴿⁱⁿᵗ[end]
    Sᴸ⁰ = zeros(N, N); Sᴿ⁰ = zeros(N, N)
    for i=1:N, j = 1:N
        Sᴸ⁰[i, j] = Sᴸⁱⁿᵗ⁻ᶠ[i, j](Rᴸ⁰);
        Sᴿ⁰[i, j] = Sᴿⁱⁿᵗ⁻ᶠ[i, j](Rᴿ⁰);
    end
    Sᴸ⁰ = round(Sᴸ⁰, 0); Sᴿ⁰ = round(Sᴿ⁰, 0)
    #@info "Target matrix values:\nSᴸ⁰=\n$(mat2string(Sᴸ⁰))\nSᴿ⁰=\n$(mat2string(Sᴿ⁰))"
    Sᴸ = Array{Function, 2}(N, N); Sᴿ = Array{Function, 2}(N, N)
    for i=1:N, j = 1:N
        if i ∈ s¹:sᵉ && j ∈ s¹:sᵉ
            Sᴸ[i, j] = R -> begin
                if R < Rᴸ⁰
                    return round(Sᴸⁱⁿᵗ⁻ᶠ[i, j](Rᴸ⁰), 0)
                elseif Rᴸ⁰ <= R < Rᴸᵉ
                    return Sᴸⁱⁿᵗ⁻ᶠ[i, j](R)
                elseif R >= Rᴸᵉ
                    return Sᴸⁱⁿᵗ⁻ᶠ[i, j](Rᴸᵉ)
                end
            end
            Sᴿ[i, j] = R -> begin
                if R < Rᴿ⁰
                    return Sᴿⁱⁿᵗ⁻ᶠ[i, j](Rᴿ⁰)
                elseif Rᴿ⁰ <= R < Rᴿᵉ
                    return Sᴿⁱⁿᵗ⁻ᶠ[i, j](R)
                elseif R >= Rᴿᵉ
                    return round(Sᴿⁱⁿᵗ⁻ᶠ[i, j](Rᴿᵉ), 0)
                end
            end
        else
            Sᴸ[i, j] = R -> Sᴸ⁰[i, j]
            Sᴿ[i, j] = R -> Sᴿ⁰[i, j]
        end
    end
    @info "Ready valuable function matrices with sizes: $(size(Sᴸ)) and $(size(Sᴿ))"
    return Sᴸ, Sᴿ
end

function internalFunctions(
    Rᴸⁱⁿᵗ::Vector{Float64}, Rᴿⁱⁿᵗ::Vector{Float64},
    Sv::Vector{Matrix{Float64}})
    Sᴸⁱⁿᵗ = Sv[1:length(Rᴸⁱⁿᵗ)]; Sᴿⁱⁿᵗ = Sv[end-length(Rᴿⁱⁿᵗ)+1:end]
    Sᴸⁱⁿᵗ⁻ᶠ = matl2matfsl(Rᴸⁱⁿᵗ, Sᴸⁱⁿᵗ, behaviour="nearest")[1]
    Sᴿⁱⁿᵗ⁻ᶠ = matl2matfsl(Rᴿⁱⁿᵗ, Sᴿⁱⁿᵗ, behaviour="nearest")[1]
    return Sᴸⁱⁿᵗ⁻ᶠ, Sᴿⁱⁿᵗ⁻ᶠ
end

function argumentGrid(Rᴸ::Vector{Float64}, Rᴿ::Vector{Float64}, Rᵛ::Vector{Float64})
    hᴸ = minimum(abs(Rᴸ[2:end] - Rᴸ[1:end-1])); Hᴸ = maximum(abs(Rᴸ[2:end] - Rᴸ[1:end-1]))
    hᴿ = minimum(abs(Rᴿ[2:end] - Rᴿ[1:end-1])); Hᴿ = maximum(abs(Rᴿ[2:end] - Rᴿ[1:end-1]))
    ΔRᴸ = abs(Rᴸ[end] - Rᴸ[1])
    ΔRᴿ = abs(Rᴿ[end] - Rᴿ[1])

    Rᴸᵉˣᵗ = Rᴸ[1] - ΔRᴸ - Hᴸ
    Rᴿᵉˣᵗ = Rᴿ[end] + ΔRᴿ + Hᴿ
    iRᴸ = findlast(R -> R <= Rᴸᵉˣᵗ, Rᵛ); iRᴸ = iRᴸ > 0 ? iRᴸ : 1
    iRᴿ = findfirst(R -> R >= Rᴿᵉˣᵗ, Rᵛ); iRᴿ = iRᴿ > 0 ? iRᴿ : length(Rᵛ)
    Rᴸᵉˣᵗ = Rᵛ[iRᴸ]; Rᴿᵉˣᵗ = Rᵛ[iRᴿ]
    @info "Smoothing intervals: left - [$Rᴸᵉˣᵗ, |$(Rᴸ[1]),$(Rᴸ[end])], right - [$(Rᴿ[1]), $(Rᴿ[end])|, $Rᴿᵉˣᵗ]"
    @info "Extra area steps: left - hᴸ = $hᴸ, right - hᴿ = $hᴿ; Max steps: left - Hᴸ = $Hᴸ, right - Hᴿ = $Hᴿ"
    return unique(vcat(collect(Rᴸᵉˣᵗ:hᴸ:Rᴸ[1]), Rᴸ)), unique(vcat(Rᴿ, collect(Rᴿ[end]:hᴿ:Rᴿᵉˣᵗ)))
end

function solutionBoundaryPoints(
    Rv::Vector{Float64}, Sv::Vector{Matrix{Float64}},
    ϵ::Float64, ϵ⁽¹⁾::Float64; rev=false, max_Δ⁽⁰⁾ = 10, max_Δ⁽¹⁾ = 10)
    max_data_count = 100
    Δ⁽⁰⁾ = 0; Δ⁽¹⁾ = 0
    R = rev == false ? Rv : reverse(Rv)
    S = rev == false ? Sv : reverse(Sv)

    kᵉ = 1; k_range = 2 : 1 : (max_data_count > length(R) ? length(R) : max_data_count)

    @info "Boundary calculation in interval R = [$(R[1]), $(R[end])], in the $(rev ? "reversed" : "direct") order."
    @info "With epsilons: ϵ = $ϵ, ϵ⁽¹⁾ = $ϵ⁽¹⁾"
    for k ∈ k_range
        kᵉ = k; kₘ = ceil(Int, (k - 1) * (1 - 1 / golden))
        kₘ = kₘ > 1 ? kₘ : 1

        ΔS⁽⁰⁾ = abs(S[k] - S[1]); ΔS⁽⁰⁾ₘ = abs(S[kₘ] - S[1])
        ΔS⁽¹⁾ = ΔS⁽⁰⁾ / abs(R[k] - R[1]); ΔS⁽¹⁾ₘ = ΔS⁽⁰⁾ₘ / abs(R[kₘ] - R[1])
        #@info "R = $(R[k]), ΔS⁽⁰⁾ = $ΔS⁽⁰⁾, ΔS⁽¹⁾ = $ΔS⁽¹⁾\n"
        if all(s -> s > 0, ΔS⁽⁰⁾) && all(s -> s > 0, ΔS⁽¹⁾)
            ϵᵏ = filter(ϵᵏᵢ -> ϵᵏᵢ > 0, ΔS⁽⁰⁾ₘ ./ ΔS⁽⁰⁾)
            ϵ⁽¹⁾ᵏ = 1 - filter(ϵ⁽¹⁾ᵏᵢ -> ϵ⁽¹⁾ᵏᵢ > 0, ΔS⁽¹⁾ₘ ./ ΔS⁽¹⁾)
            @info "k = $k; R = $(R[k]); ϵᵏ = $ϵᵏ <-> ϵ⁽¹⁾ᵏ = $ϵ⁽¹⁾ᵏ; ΔS⁽⁰⁾ = $ΔS⁽⁰⁾; ΔS⁽⁰⁾ₘ = $ΔS⁽⁰⁾ₘ"
            if any(ϵᵏᵢ -> ϵᵏᵢ > ϵ, ϵᵏ)
                Δ⁽⁰⁾ += 1
            end
            if any(ϵ⁽¹⁾ᵏᵢ -> ϵ⁽¹⁾ᵏᵢ > ϵ⁽¹⁾, ϵ⁽¹⁾ᵏ)
                Δ⁽¹⁾ += 1
            end
        end
        if Δ⁽⁰⁾ > max_Δ⁽⁰⁾ && Δ⁽¹⁾ > max_Δ⁽¹⁾
            @warn "Count reached for k = $kᵉ at R = $(R[kᵉ])"
            break
        end
    end
    if rev == false
        @info "The result range is [$(R[1]), $(R[kᵉ])]"
        return R[1:kᵉ]
    else
        @info "The result range is [$(R[kᵉ]), $(R[1])]"
        return reverse(R[1:kᵉ])
    end
end

dS_dR₂ = (Sᵏ⁺¹, Sᵏ⁻¹, h) -> begin return (Sᵏ⁺¹ - Sᵏ⁻¹) / (2 * h) end
dS_dR₄ = (Sᵏ⁺², Sᵏ⁺¹, Sᵏ⁻¹, Sᵏ⁻², h) -> begin return (-Sᵏ⁺² + 8 * Sᵏ⁺¹ - 8 * Sᵏ⁻¹ + Sᵏ⁻²) / (12 * h) end

# Sᵏ⁻² = view(S[k - 2], s¹:sᵉ, s¹:sᵉ)
# Sᵏ⁻¹ = view(S[k - 1], s¹:sᵉ, s¹:sᵉ)
# Sᵏ   = view(S[k], s¹:sᵉ, s¹:sᵉ)
# Sᵏ⁺¹ = view(S[k + 1], s¹:sᵉ, s¹:sᵉ)
# Sᵏ⁺² = view(S[k + 2], s¹:sᵉ, s¹:sᵉ)
#
# Rᵏ = k > 2 ? R[k-2:k+2] : R[k-1:k+1]
# h = mean(Rᵏ[2:end] - Rᵏ[1:end-1]); push!(steps, h) # rough step size for good functions
# h = rev == false ? h : -h
#
# S⁽¹⁾ᵏ = k > 2 ? dS_dR₄(Sᵏ⁺², Sᵏ⁺¹, Sᵏ⁻¹, Sᵏ⁻², h) : dS_dR₂(Sᵏ⁺¹, Sᵏ⁻¹, h)
#
#
# ΔS⁽⁰⁾ᵏ = Sᵏ - Sᵏ⁻¹
# ΔS⁽¹⁾ᵏ = k > 2 ? S⁽¹⁾ᵏ - S⁽¹⁾ᵏ⁻¹ : zeros(N, N)
# push!(vΔS⁽⁰⁾, ΔS⁽⁰⁾ᵏ); push!(vΔS⁽¹⁾, ΔS⁽¹⁾ᵏ)

# S⁽¹⁾ᵏ⁻¹ = zeros(N, N); S⁽¹⁾ᵏ = zeros(N, N)
# vΔS⁽⁰⁾ = Vector{Float64}(undef, 0); ΔS⁽⁰⁾ᵏ = zeros(N, N)
# vΔS⁽¹⁾ = Vector{Float64}(undef, 0); ΔS⁽¹⁾ᵏ = zeros(N, N)


function calculate∂²_∂R²(Rᵖᵒⁱⁿᵗˢ::Vector{Float64},
    ∂_∂Rᴰᵈᵃᵗᵃ::Matrix{Float64},
    N::Int)
    @assert length(Rᵖᵒⁱⁿᵗˢ) == size(∂_∂Rᴰᵈᵃᵗᵃ, 1)

    L = length(Rᵖᵒⁱⁿᵗˢ)
    M = size(∂_∂Rᴰᵈᵃᵗᵃ, 2)
    @assert M == dataSizeOfSymetricMatrix(N) "$M≠$(dataSizeOfSymetricMatrix(N))"

    ∂²_∂R²ᴰᵈᵃᵗᵃ = Matrix{Float64}(undef, L, M)
    ∂²_∂R²ᴰᵈᵃᵗᵃ_diag = Matrix{Float64}(undef, L, N)
    ∂_∂Rᴰᵈᵃᵗᵃ_func = Array{Function, 2}(undef, N, N)
    ∂_∂Rᴰᵈᵃᵗᵃ_spl = Array{Dierckx.Spline1D, 2}(undef, N, N)
    lᵖ = 1
    for i=1:N, j=i+1:N
      l = dataColumnOfSymetricMatrix(i, j, N)
      @assert lᵖ <= l "$lᵖ>$l"
      X = Rᵖᵒⁱⁿᵗˢ; Y = ∂_∂Rᴰᵈᵃᵗᵃ[:, l]
      spl = Dierckx.Spline1D(X, Y; w=ones(length(X)), k=2, bc="nearest", s=0.0)
      ispl = Dierckx.Spline1D(X, -Y; w=ones(length(X)), k=2, bc="nearest", s=0.0)

      ∂_∂Rᴰᵈᵃᵗᵃ_func[i, j] = R -> Dierckx.evaluate(spl, R)
      ∂_∂Rᴰᵈᵃᵗᵃ_func[j, i] = R -> Dierckx.evaluate(ispl, R)

      ∂_∂Rᴰᵈᵃᵗᵃ_spl[i, j] = spl; ∂_∂Rᴰᵈᵃᵗᵃ_spl[j, i] = ispl
    end
    for i=1:N
      ∂_∂Rᴰᵈᵃᵗᵃ_func[i, i] = R -> 0.0
      ∂_∂Rᴰᵈᵃᵗᵃ_spl[i, i] = Dierckx.Spline1D(Rᵖᵒⁱⁿᵗˢ, zeros(L); w=ones(length(Rᵖᵒⁱⁿᵗˢ)), k=1, bc="nearest", s=0.0)
    end

    ∂²_∂R²ᴰᵈᵃᵗᵃ = Matrix{Float64}(undef, L, M)
    ∂²_∂R²ᴰᵈᵃᵗᵃ_diag = Matrix{Float64}(undef, L, N)
    for lʳ = 1:L
      R = Rᵖᵒⁱⁿᵗˢ[lʳ]
      τ⁽¹⁾ = matf2mat(R, ∂_∂Rᴰᵈᵃᵗᵃ_func)
      ∇τ⁽¹⁾ = Dierckx.derivative.(∂_∂Rᴰᵈᵃᵗᵃ_spl, R; nu=1)
      τ⁽²⁾ = τ⁽¹⁾*τ⁽¹⁾ + ∇τ⁽¹⁾ # ???

      ∂²_∂R²ᴰᵈᵃᵗᵃ[lʳ, :] = matl2matupper(∇τ⁽¹⁾)
      ∂²_∂R²ᴰᵈᵃᵗᵃ_diag[lʳ, :] = diag(τ⁽²⁾)
    end

    return ∂²_∂R²ᴰᵈᵃᵗᵃ, ∂²_∂R²ᴰᵈᵃᵗᵃ_diag
end
