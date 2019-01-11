import Base
import Base.-
using NumericalIntegration
using DataStructures

mutable struct OptimizableArea
    ∂_∂R::Function
    area::LandauZenerArea
    limits::Tuple{Float64, Float64}
    X_data::Vector{Float64}
    Y_data::Vector{Float64}
    ∫∂_∂RdR::Float64
    ϵ::Float64
end
function -(oa::OptimizableArea)
    return OptimizableArea(R -> -oa.∂_∂R(R), -oa.area, oa.limits, oa.X_data, -oa.Y_data, -oa.∫∂_∂RdR, oa.ϵ)
end
function Base.show(io::IO, mime::AbstractString, oa::OptimizableArea)
    show(io, oa)
end
function Base.show(io::IO, oa::OptimizableArea)
    i = oa.area.states[1]; j = oa.area.states[2]
    print(io, "Ω{⦗$(i)|∂/∂R|$(j)⦘, ∫⟨$(i)|∂/∂R|$(j)⟩dR = $(oa.∫∂_∂RdR) in [$(oa.limits[1]), $(oa.limits[2])], ϵ = $(oa.ϵ)}")
end
function Base.show(io::IO, mime::AbstractString, voa::Vector{OptimizableArea})
    show(io, voa)
end
function Base.show(io::IO, voa::Vector{OptimizableArea})
    print(io, "⦃$(join(voa, "; "))⦄")
end
function Base.show(io::IO, mime::AbstractString, moa::Matrix{Vector{OptimizableArea}})
    show(io, moa)
end
function Base.show(io::IO, moa::Matrix{Vector{OptimizableArea}})
    println(io, "Optimizable areas:")
    for i = 1:size(moa, 1), j = 1:size(moa, 2)
        if isassigned(moa, i, j) && !isempty(moa[i, j])
            println(io, "  ⦿ ⟨$(i)|∂/∂R|$(j)⟩: $(moa[i, j])")
        end
    end
end

function getOptimizedCouplings(
    oam::Matrix{Vector{OptimizableArea}},
    ∂_∂R::Matrix{Function}, ∂_∂Rᵐᵒᵈᵉˡ::Matrix{Function})
    @assert size(∂_∂R)[1] == size(∂_∂R)[2]
    N = size(∂_∂R)[1]
    ∂_∂Rᴼ = Matrix{Function}(undef, N, N)
end

function findOptimizableCouplings(
    ∂_∂R::Matrix{Function}, ∂_∂R_data::Matrix{Float64},
    areas::Matrix{Vector{LandauZenerArea}},
    solutions::Vector{LocalSolution},
    nonadiabatic_config::NonadiabaticAreasConfiguration)
    @info "⟑⟑⟑Searching for a better optimizable couplings... in matrix $(size(∂_∂R)[1])×$(size(∂_∂R)[2])"
    @assert size(∂_∂R) == size(areas)
    @assert size(∂_∂R)[1] == size(∂_∂R)[2]
    @assert dataSizeOfSymetricMatrix(size(∂_∂R)[1]) == size(∂_∂R_data)[2] - 1

    Α_config = nonadiabatic_config.nonadiabatic_areas[SINGLE_PEAK::NonadiabaticAreaTypes]
    ϵˡᶻ = abs(Α_config.error_potential_∂_∂R_coordinate)

    N = size(∂_∂R)[1]
    R_data = ∂_∂R_data[:, 1]
    optim_areas = Matrix{Vector{OptimizableArea}}(undef, N, N)
    for i ∈ 1:N, j ∈ 1:N
        if i < j
            if !isassigned(areas, i, j) || isempty(areas[i, j])
                @info "No data for ⟨$(i)|∂/∂R|$(j)⟩... skipping"
                continue
            end
            all_solution_states = states_of_solution((i, j), solutions)
            if length(all_solution_states) > 2
                @warn "Skipping ⟨$(i)|∂/∂R|$(j)⟩ with solution states $(all_solution_states) > 2"
                continue
            end
            @info "Peaks of ⟨$(i)|∂/∂R|$(j)⟩..."
            optim_areas_ij = Vector{OptimizableArea}()
            for k ∈ 1:length(areas[i, j])
                area = areas[i, j][k]
                R₀ = area.R₀
                indexes_R₀ = findall(R -> abs(R - R₀) ≤ ϵˡᶻ, R_data)
                @assert !isempty(indexes_R₀)

                l = dataColumnOfSymetricMatrix(i, j, N)
                ∂_∂R_ij = ∂_∂R_data[:, l + 1]

                ixR₀ = find_closest_index(R₀, R_data, ϵˡᶻ)
                iter_step = give_iterations_and_step(R_data, area.Rₐ, area.Rᵦ)
                iterations = round(Int, 2 * golden * iter_step[1]); ΔR = iter_step[2]
                iterations = iterations > 161 ? iterations : 161
                ∫∂_∂R_ij = 0.0; ∫∂_∂R_model = π/2;  ϵ_∫ = 1e-2
                @info "Looking at the ∫⟨$(i)|∂/∂R|$(j)⟩dR in the area [$(area.Rₐ), $(area.Rᵦ)](ixR₀ = $(ixR₀)) with ΔR = $(ΔR) and $iterations iterations..."
                ΔI = 1
                ix_l = ixR₀-ΔI ≥ 1 ? ixR₀-ΔI : 1;
                ix_r = ixR₀+ΔI ≤ length(R_data) ? ixR₀+ΔI : length(R_data);
                for n_iter ∈ 1:iterations
                    ΔI = n_iter
                    ix_l = ixR₀-ΔI ≥ 1 ? ixR₀-ΔI : 1;
                    ix_r = ixR₀+ΔI ≤ length(R_data) ? ixR₀+ΔI : length(R_data);
                    X_R_data = R_data[ix_l:ix_r]
                    Y_∂_∂R = ∂_∂R_ij[ix_l:ix_r]
                    @assert length(X_R_data) == length(Y_∂_∂R)
                    ∫∂_∂R_ij = NumericalIntegration.integrate(X_R_data, Y_∂_∂R)
                    if abs(abs(∫∂_∂R_ij) - ∫∂_∂R_model) ≤ ϵ_∫ || abs(∫∂_∂R_ij) ≥ ∫∂_∂R_model
                        @info "|∫⟨$(i)|∂/∂R|$(j)⟩dR| = $(abs(∫∂_∂R_ij)) ≈ π/2 in the interval [$(X_R_data[1]), $(X_R_data[end])], δ=$(4*area.τ₀), ⤒ = $(1/(4*area.τ₀)); iteration $n_iter of $iterations."
                        optim_area = OptimizableArea(∂_∂R[i, j], area, (X_R_data[1], X_R_data[end]), X_R_data, Y_∂_∂R, ∫∂_∂R_ij, ϵ_∫)
                        push!(optim_areas_ij, optim_area)
                        break
                    end
                end
                if abs(abs(∫∂_∂R_ij) - ∫∂_∂R_model) > ϵ_∫
                    @warn "|∫⟨$(i)|∂/∂R|$(j)⟩dR| = $(abs(∫∂_∂R_ij)) ≠ $(π/2) ≠ π/2 in the interval [$(R_data[ix_l]), $(R_data[ix_r])], δ=$(4*area.τ₀), ⤒ = $(1/(4*area.τ₀)); with $iterations iterations."
                end
                @info "... Done."
            end
            optim_areas[i, j] = optim_areas_ij
            optim_areas[j, i] = map(oa -> -oa, optim_areas_ij)
            # ----
        end
    end
    return optim_areas
end

function optimizeCouplings(
    oam::Matrix{Vector{OptimizableArea}},
    ∂_∂R::Matrix{Function}, ∂_∂R_model::Matrix{Function})
    @assert size(∂_∂R) == size(∂_∂R_model)
    @assert size(oam) == size(∂_∂R)
    @assert size(∂_∂R, 1) == size(∂_∂R, 2)
    N = size(∂_∂R, 1)
    ∂_∂R_model_optimized = Matrix{Function}(undef, N, N)
    for i = 1:N, j = 1:N
        if i == j
            ∂_∂R_model_optimized[i, j] = R -> 0.0
        elseif i < j
            if isassigned(oam, i, j) && !isempty(oam[i, j])
                @info "Making optimized model coupling ⟨$(i)|∂/∂R|$(j)⟩..."
                oav = oam[i, j]
                definition_list = Vector{Pair{Tuple{Float64, Float64}, Tuple{Function, Bool}}}(undef, 0)
                original_functionⁱʲ = ∂_∂R[i, j]
                model_functionⁱʲ = ∂_∂R_model[i, j]
                @info "Function intervals defining..."
                for k = 1:length(oav)
                    oaₖⁱʲ = oav[k]
                    limitsₖⁱʲ = oaₖⁱʲ.limits
                    if k == 1
                        definition₋₁ = (-Inf, limitsₖⁱʲ[1]) => (model_functionⁱʲ, false)
                        push!(definition_list, definition₋₁)
                        @info "⇒𝔹 Function interval (-∞, $(limitsₖⁱʲ[1])) → ⟨$(i)|∂/∂R|$(j)⟩ᴹᴼᴰᴱᴸ"
                    elseif k > 1
                        oaₖ₋₁ⁱʲ = oav[k - 1]
                        limitsₖ₋₁ⁱʲ = oaₖ₋₁ⁱʲ.limits
                        definition₋₁ = (limitsₖ₋₁ⁱʲ[2], limitsₖⁱʲ[1]) => (model_functionⁱʲ, false)
                        push!(definition_list, definition₋₁)
                        @info " - Function interval ($(limitsₖ₋₁ⁱʲ[2]), $(limitsₖⁱʲ[1])) → ⟨$(i)|∂/∂R|$(j)⟩ᴹᴼᴰᴱᴸ"
                    else
                        # nothing
                    end
                    definition = (limitsₖⁱʲ[1], limitsₖⁱʲ[2]) => (original_functionⁱʲ, true)
                    push!(definition_list, definition)
                    @info " - Function interval ($(limitsₖⁱʲ[1]), $(limitsₖⁱʲ[2])) → ⟨$(i)|∂/∂R|$(j)⟩"
                    if k == length(oav)
                        definition₊₁ = (limitsₖⁱʲ[2], Inf) => (model_functionⁱʲ, false)
                        push!(definition_list, definition₊₁)
                        @info "⇒𝔼 Function interval ($(limitsₖⁱʲ[2]), ∞) → ⟨$(i)|∂/∂R|$(j)⟩ᴹᴼᴰᴱᴸ"
                    end
                end
                definitions = OrderedDict{Tuple{Float64, Float64}, Tuple{Function, Bool}}(definition_list)
                @info "Piecewise function defining with definitions:\n$definitions"
                fⁱʲ, pⁱʲ = piecewise_function(definitions)
                ∂_∂R_model_optimized[i, j] = fⁱʲ
                ∂_∂R_model_optimized[j, i] = R-> -fⁱʲ(R)
                @info "Breakpoints are defined: $pⁱʲ"
                @info "... Done for ⟨$(i)|∂/∂R|$(j)⟩."
            else
                ∂_∂R_model_optimized[i, j] = ∂_∂R_model[i, j]
                ∂_∂R_model_optimized[j, i] = ∂_∂R_model[j, i]
            end
        end
    end
    @info "... Done."
    return ∂_∂R_model_optimized
end


function states_of_solution(braket::Tuple{Int, Int}, solutions::Vector{LocalSolution})
    states = Vector{Int}(undef, 0)
    sort(unique(convert(Vector{Int},
        flatten(map(s -> s.states,
            filter(s -> braket[1] ∈ s.states && braket[2] ∈ s.states, solutions))))))
end

function give_iterations_and_step(R_data::Vector{Float64}, Rᵃ::Float64, Rᵇ::Float64)
    @assert Rᵃ < Rᵇ
    R_data_interval = filter(R -> Rᵃ ≤ R ≤ Rᵇ, R_data)
    R_data_right = R_data_interval[2:end]
    R_data_left = R_data_interval[1:end-1]
    @assert length(R_data_left) == length(R_data_right)
    set_ΔR = abs.(R_data_right - R_data_left)
    ΔR = minimum(set_ΔR)
    @assert ΔR ≠ 0
    return (round(Int, (Rᵇ - Rᵃ) / ΔR), ΔR)
end

function find_closest_index(R₀::Float64, R_data::Vector{Float64}, ϵᴿ::Float64)
    @info "Find closest index with ϵᴿ = $ϵᴿ for R₀ = $R₀"
    pR = mapreduce(pR -> (pR[1], abs(pR[2] - R₀)),
            (p1, p2) -> p1[2] ≤ p2[2] ? p1 : p2,
            collect(enumerate(R_data)))
    R = R_data[pR[1]]
    ΔR = abs(R - R₀)
    @info "Found minimum distance point R = $R with ΔR = $ΔR"
    return pR[1]
end
