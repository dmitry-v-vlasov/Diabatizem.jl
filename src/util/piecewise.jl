using DataStructures

const 𝐆old = MathConstants.golden

function piecewise_function(
  definitions::OrderedDict{Tuple{Float64, Float64}, Tuple{Function, Bool}}; ϵᵇᵖ::Float64=1e-2)
  @assert all(d -> d.first[1] < d.first[2], definitions)
  keypoints = flatten(keys(definitions))
  @assert issorted(keypoints)
  inf_point_pos = sort(findall(keypoint -> abs(keypoint) == Inf, keypoints))
  @assert length(inf_point_pos) ≤ 2
  @assert all(pos -> pos ∈ [1, length(keypoints)], inf_point_pos)
  @assert ifelse(length(inf_point_pos) == 2, keypoints[1] == -Inf, true)
  @assert ifelse(length(inf_point_pos) == 2, keypoints[end] == Inf, true)
  @assert ifelse(length(inf_point_pos) == 1 && inf_point_pos[1] == 1, keypoints[1] == -Inf, true)
  @assert ifelse(length(inf_point_pos) == 1 && inf_point_pos[1] == length(keypoints), keypoints[end] == Inf, true)
  @assert ifelse(length(inf_point_pos) == 0, abs(keypoints[1]) ≠ Inf, true)
  @assert ifelse(length(inf_point_pos) == 0, abs(keypoints[end]) ≠ Inf, true)

  lfuncs = map(d -> d.first[2] => d.second, collect(definitions)[1:end-1])
  rfuncs = map(d -> d.first[1] => d.second, collect(definitions)[2:end])
  @assert length(lfuncs) == length(rfuncs)
  Nₚ = length(lfuncs)
  point_pair_funcs = Vector{Pair{Tuple{Float64,Float64}, Tuple{Tuple{Function, Bool},Tuple{Function, Bool}}}}()
  for k = 1:Nₚ
    ppf = (lfuncs[k].first, rfuncs[k].first) => (lfuncs[k].second, rfuncs[k].second)
    push!(point_pair_funcs, ppf)
  end
  @assert all(ppf -> abs(ppf.first[1]) ≠ Inf, point_pair_funcs)
  @assert all(ppf -> abs(ppf.first[2]) ≠ Inf, point_pair_funcs)
  @assert all(ppf -> abs(ppf.first[1] - ppf.first[2]) ≤ ϵᵇᵖ, point_pair_funcs)
  breakpoint_funcs = map(ppf -> (ppf.first[1] + ppf.first[2]) / 2 => ppf.second, point_pair_funcs)
  breakpoint_func_gaps = map(bf -> bf.first => (bf.second[2][1](bf.first + ϵᵇᵖ) - bf.second[1][1](bf.first - ϵᵇᵖ), bf.second[1][2], bf.second[2][2]), breakpoint_funcs)
  @assert length(breakpoint_funcs) == length(breakpoint_func_gaps)
  breakpoint_parameters = map(
  bp_entry -> begin
      k = bp_entry[1]
      func_pars = bp_entry[2]
      func_gap = func_pars.second[1]
      keep_left = func_pars.second[2]
      keep_right = func_pars.second[3]
      Δy = func_gap; Δx = abs(𝐆old * Δy)
      Δx = Δx ≥ ϵᵇᵖ ? Δx : ϵᵇᵖ
      x₀ = func_pars.first
      f1 = breakpoint_funcs[k].second[1][1]
      f2 = breakpoint_funcs[k].second[2][1]
      # @assert typeof(f1) == Function
      # @assert typeof(f2) == Function
      return (x₀, Δx, Δy, (f1, f2), keep_left, keep_right)
    end,
    enumerate(breakpoint_func_gaps))
  return make_piecewise_sigmoid_function(definitions, breakpoint_funcs, breakpoint_parameters)
end

function make_piecewise_sigmoid_function(
  definitions::OrderedDict{Tuple{Float64, Float64},Tuple{Function, Bool}},
  breakpoint_funcs, breakpoint_parameters;
  ϵᵇᵖ::Float64=1e-2)
  @assert length(definitions) ≥ 2

  breakpoints = Vector{Float64}(undef, 0)
  functions = Vector{Function}(undef, 0)
  definition_list = collect(definitions)
  for k = 1:length(definition_list)-1
    defₖ_left = definition_list[k]
    xₖ_a_left = defₖ_left.first[1]; xₖ_b_left = defₖ_left.first[2]
    @assert xₖ_a_left < xₖ_b_left
    @assert ifelse(xₖ_a_left == -Inf, k == 1, true)
    keep_leftₖ = defₖ_left.second[2]
    Δxₖ_left = xₖ_b_left - xₖ_a_left
    if keep_leftₖ
      @warn "Interval [$xₖ_a_left, $xₖ_b_left] will be kept (untouched function)."
    end

    defₖ_right = definition_list[k + 1]
    xₖ_a_right = defₖ_right.first[1]; xₖ_b_right = defₖ_right.first[2]
    @assert xₖ_a_right < xₖ_b_right
    @assert ifelse(xₖ_b_right == Inf, k + 1 == length(definition_list), true)
    keep_rightₖ = defₖ_right.second[2]
    Δxₖ_right = xₖ_b_right - xₖ_a_right
    if keep_rightₖ
      @warn "Interval [$xₖ_a_right, $xₖ_b_right] will be kept (untouched function)."
    end

    parametersₖ = breakpoint_parameters[k]
    x₀ₖ = parametersₖ[1]
    Δxₛ = parametersₖ[2]
    Δyₛ = parametersₖ[3]
    f1, f2 = parametersₖ[4]
    kf1 = parametersₖ[5]
    kf2 = parametersₖ[6]
    @assert kf1 == keep_leftₖ
    @assert kf2 == keep_rightₖ

    δxₖ_left = ifelse(xₖ_a_left == -Inf,
      10 * 𝐆old * Δxₛ,
      ifelse((𝐆old - 1) * Δxₛ < Δxₖ_left, (𝐆old - 1) * Δxₛ, (2 - 𝐆old) * Δxₖ_left))

    δxₖ_right = ifelse(xₖ_b_right == Inf,
      10 * 𝐆old * Δxₛ,
      ifelse((𝐆old - 1) * Δxₛ < Δxₖ_right, (𝐆old - 1) * Δxₛ, (2 - 𝐆old) * Δxₖ_right))
    #@assert(δxₖ_left ≥ (𝐆old - 1) * Δxₛ / 2, "δxₖ_left ≥ (𝐆old - 1)Δxₛ/2: δxₖ_left = $(δxₖ_left), (𝐆old - 1)Δxₛ/2 = $((𝐆old - 1)*Δxₛ/2)")
    #@assert(δxₖ_right ≥ (𝐆old - 1) * Δxₛ / 2, "δxₖ_right ≥ (𝐆old - 1)Δxₛ/2: δxₖ_right = $(δxₖ_right),  (𝐆old - 1)Δxₛ/2 = $((𝐆old - 1)*Δxₛ/2)")

    α = 0.2 * (δxₖ_left + δxₖ_right)
    δxₖ_left_new = keep_leftₖ ? δxₖ_left / 10.0 : δxₖ_left
    δxₖ_right_new = keep_rightₖ ? δxₖ_right / 10.0 : δxₖ_right
    if δxₖ_left_new < δxₖ_left
      x₀ₖ = x₀ₖ + δxₖ_left - δxₖ_left_new
    end
    if δxₖ_right_new < δxₖ_right
      x₀ₖ = x₀ₖ - δxₖ_right + δxₖ_right_new
    end
    @info "Blending function $f1; $f2; $x₀ₖ; $α"
    blending_functionₖ = sigmoid_of_name(f1, f2, x₀ₖ, α, "")

    breakpoint_left = x₀ₖ - δxₖ_left_new
    breakpoint_right = x₀ₖ + δxₖ_right_new
    function_left = defₖ_left.second[1]
    function_middle = blending_functionₖ
    function_right = defₖ_right.second[1]

    @info "Made blending function in interval:
    [$(breakpoint_left), $(breakpoint_right)],
    α = $α, x₀ₖ = $x₀ₖ, δxₖ_left = $δxₖ_left_new, δxₖ_right = $δxₖ_right_new,
    keep_leftₖ = $keep_leftₖ, keep_rightₖ = $keep_rightₖ,
    δxₖ_left_old = $δxₖ_left, δxₖ_right_old = $δxₖ_right"

    push!(functions, function_left)
    push!(breakpoints, breakpoint_left)
    push!(functions, function_middle)
    push!(breakpoints, breakpoint_right)
    if k + 1 == length(definition_list)
      push!(functions, function_right)
    end
  end
  pwf = make_piecewise_function(breakpoints, functions)
  return pwf, breakpoints
end

function make_piecewise_function(
  breakpoints::Vector{Float64}, functions::Vector{Function})
  pwf(R) = begin
    nᶠ = searchsortedfirst(breakpoints, R)
    f_piecewise = functions[nᶠ]
    return f_piecewise(R)
  end
end

const sigmoid = (x, x₀, α) -> 1/(1 + exp(-(x - x₀)/α))
function sigmoid_of_name(f1::Function, f2::Function, x₀, α, name::AbstractString)
  return x -> begin
    sf = (1 - sigmoid(x, x₀, α))*f1(x) + sigmoid(x, x₀, α)*f2(x)
    if sf === NaN || sf === NaN64 || sf === NaN32 || sf === NaN16
      @error "SIGMOID NaN: x=$x, x₀=$x₀, α=$α, f1(x)=$(f1(x)), f2(x)=$(f2(x)), σ(x)=$(sigmoid(x, x₀, α)), name=$name"
    end
    return sf
  end
end
