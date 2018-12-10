using JSON
using DataFrames
using DataStructures
using Nullables
using CSVFiles

mutable struct InputPaths
  file_hamiltonian_adiabatic::AbstractString
  file_coupling_∂_∂R_adiabatic::AbstractString
  file_coupling_∂_∂R_adiabatic_model::Nullable{AbstractString}
  file_transformation_matrix_initial::Nullable{AbstractString}
  file_transformation_matrix::Nullable{AbstractString}
end

mutable struct InputData
  hamiltonian_adiabatic::DataFrame
  coupling_∂_∂R_adiabatic::DataFrame
  coupling_∂_∂R_adiabatic_model::Nullable{DataFrame}
end

mutable struct OutputPaths
  file_potentials_diabatic::AbstractString
  file_hamiltonian_diabatic::AbstractString
  file_coupling_∂_∂R_diabatic::AbstractString
  file_coupling_∂²_∂R²_diabatic::AbstractString
  file_coupling_∂²_∂R²_diabatic_diag::AbstractString
  file_transformation_matrix::AbstractString
  file_transformation_matrix_general::AbstractString
end

@enum CalculationStrategy LANDAU_ZENER LANDAU_ZENER_WITH_EXTERNAL_MODEL_DATA EXTERNAL_MODEL_DATA JUST_MAKE_THEM_ZEROS UNKNOWN
const MAPPING_CalculationStrategy = Dict(
  "Landau-Zener" => LANDAU_ZENER::CalculationStrategy,
  "Landau-Zener-with-external-model-data" => LANDAU_ZENER_WITH_EXTERNAL_MODEL_DATA::CalculationStrategy,
  "external-model-data" => EXTERNAL_MODEL_DATA::CalculationStrategy,
  "just-make-them-zeros" => JUST_MAKE_THEM_ZEROS::CalculationStrategy
)

@enum NonadiabaticAreaTypes SINGLE_PEAK DOUBLE_PEAK UNSUPPORTED
const MAPPING_NonadiabaticAreaTypes = Dict(
  "single-peak" => SINGLE_PEAK::NonadiabaticAreaTypes
  #"double-peak" => DOUBLE_PEAK
)

mutable struct UtilitySettings
  channel_ionic_number::Int
  channel_lowest_number::Int
end

@enum InterpolationType SPLINE_LINEAR SPLINE_QUADRATIC SPLINE_CUBIC UNSUPPORTED_INTERP
const MAPPING_InterpolationType = Dict(
  "spline-linear" => SPLINE_LINEAR::InterpolationType,
  "spline-quadratic" => SPLINE_QUADRATIC::InterpolationType,
  "spline-cubic" => SPLINE_CUBIC::InterpolationType
)
mutable struct InterpolationSettings
  hamiltonian::InterpolationType
  coupling_∂_∂R::InterpolationType
end

mutable struct AsymptoticSettings
  coordinate_start::Float64
  coordinate_step::Float64
  coordinate_step_error::Float64
  coordinate_safety_step::Float64
  potential_asymptotic_value_error::Float64
  ∂_∂R_asymptotic_value_error::Float64
  ∂_∂R_zero_value_error::Float64
end

abstract type NonadiabaticAreaSettings end

mutable struct SinglePeakNonadiabaticAreaSettings <: NonadiabaticAreaSettings
  error_∂_∂R_peak::Float64
  vanishing_∂_∂R_value::Float64
  error_vanishing_∂_∂R_value::Float64
  error_potential_distance_minimal::Float64
  error_potential_distance_coordinate::Float64
  error_potential_∂_∂R_coordinate::Float64
end
mutable struct NonadiabaticAreasConfiguration
  coordinate_start::Float64
  coordinate_step::Float64
  coordinate_piece::Float64
  coordinate_step_error::Float64
  nonadiabatic_areas::Dict{NonadiabaticAreaTypes, NonadiabaticAreaSettings}
end

mutable struct SelectedDiabatizationArea
    coordinate::Float64
    states::Tuple{Int, Int}
    extra_length::Tuple{Float64, Float64}
    bunch_exclude::Bool
end

mutable struct DiabatizationSettings
  areas::Vector{SelectedDiabatizationArea}
  keep_initial_conditions::Bool
  coordinate_start::Float64
  coordinate_step::Tuple{Float64, Float64}
  coordinate_stop::Float64
  area_closeness::Float64
  use_last_transformation_matrix_from::Nullable{Float64}
end

mutable struct CalculationSettings
  strategy::CalculationStrategy
  coordinate_start::Float64
  coordinate_step::Float64
  coordinate_piece::Float64
  coordinate_step_error::Float64
  asymptotics::AsymptoticSettings
  nonadiabatic_areas::NonadiabaticAreasConfiguration
  diabatization::DiabatizationSettings
  utility::UtilitySettings
  interpolation::InterpolationSettings
end

mutable struct Configuration
  input_paths::InputPaths
  input_data::InputData
  output_paths::OutputPaths
  settings::CalculationSettings
end

function loadConfiguration(filePath::AbstractString)
  js = JSON.parsefile(filePath; dicttype=DataStructures.OrderedDict)

  input_paths = InputPaths(
    js["input-data"]["hamiltonian-adiabatic"],
    js["input-data"]["coupling_∂_∂R_adiabatic"],
    Nullable{AbstractString}(),
    haskey(js["input-data"], "transformation-matrix-initial") ?
      js["input-data"]["transformation-matrix-initial"] :
      Nullable{AbstractString}(),
    haskey(js["input-data"], "transformation-matrix") ?
      js["input-data"]["transformation-matrix"] :
      Nullable{AbstractString}())

  Hₐ_data = loadRawData(input_paths.file_hamiltonian_adiabatic)
  fixPotentialAsymptotics!(Hₐ_data, loadUtilitySettings(js["settings"]["utility"]))
  ∂_∂Rᴬ_data = loadRawData(input_paths.file_coupling_∂_∂R_adiabatic)
  input_data = InputData(Hₐ_data, ∂_∂Rᴬ_data, Nullable{DataFrame}())

  output_paths = OutputPaths(
    js["output-data"]["potentials-diabatic"],
    js["output-data"]["hamiltonian-diabatic"],
    js["output-data"]["coupling-∂_∂R-diabatic"],
    js["output-data"]["coupling-∂²_∂R²-diabatic"],
    js["output-data"]["coupling-∂²_∂R²-diabatic-diagonal"],
    js["output-data"]["transformation-matrix"],
    js["output-data"]["transformation-matrix-general"])

  settings = loadCalculationSettings(js, input_paths, input_data)

  return Configuration(input_paths, input_data, output_paths, settings)
end

function fixPotentialAsymptotics!(Hₐ_data::DataFrame, utilitySettings)
  if utilitySettings.channel_lowest_number < 0
    @info "Fixing channel asymptotics has been skept in the configuration loading stage."
    return
  end
  Nₚ = size(Hₐ_data, 1)
  Nᵩ = size(Hₐ_data, 2) - 1
  lowest_channel = utilitySettings.channel_lowest_number > 0 ? utilitySettings.channel_lowest_number : 1
  Uₗ_∞ = Hₐ_data[Nₚ, lowest_channel + 1]
  @info "Fixing chennel asymptotics with V₀ᴬ(∞)=$(Uₗ_∞)."
  for line = 1:Nₚ, channel = 1:Nᵩ
    Hₐ_data[line, channel + 1] = Hₐ_data[line, channel + 1] - Uₗ_∞
  end
end

function loadCalculationSettings(js, input_paths::InputPaths, input_data::InputData)
  jss = js["settings"]
  strategy_name = CalculationStrategy(deriveStrategy(jss["strategy"]))
  if LANDAU_ZENER_WITH_EXTERNAL_MODEL_DATA::CalculationStrategy == strategy_name || EXTERNAL_MODEL_DATA::CalculationStrategy == strategy_name
    if haskey(js["input-data"], "coupling_∂_∂R_adiabatic_model")
      input_paths.file_coupling_∂_∂R_adiabatic_model = js["input-data"]["coupling_∂_∂R_adiabatic_model"]
    else
      throw(DomainError("The configuration setting 'coupling_∂_∂R_adiabatic_model' is obligatory for the strategy $strategy_name"))
    end
  end
  if LANDAU_ZENER_WITH_EXTERNAL_MODEL_DATA::CalculationStrategy == strategy_name || EXTERNAL_MODEL_DATA::CalculationStrategy == strategy_name
    input_data.coupling_∂_∂R_adiabatic_model = loadRawData(get(input_paths.file_coupling_∂_∂R_adiabatic_model))
  end

  coordinate_start = haskey(jss, "coordinate-start") ? (jss["coordinate-start"] < 0 ? -1 : jss["coordinate-start"]) : -1
  coordinate_step = haskey(jss, "coordinate-step") ? (jss["coordinate-step"] < 0 ? -1 : jss["coordinate-step"]) : -1
  coordinate_piece = haskey(jss, "coordinate-piece") ? (jss["coordinate-piece"] < 0 ? -1 : jss["coordinate-piece"]) : -1
  coordinate_step_error = haskey(jss, "coordinate-step-error") ? (jss["coordinate-step-error"] < 0 ? -1 : jss["coordinate-step-error"]) : -1

  asymptoticSettings = loadAsymptotics(jss["asymptotics"], coordinate_start, coordinate_step, coordinate_step_error)
  nonadiabaticAreas = loadNonadiabaticAreas(jss["nonadiabatic-areas"], coordinate_start, coordinate_step, coordinate_piece, coordinate_step_error)
  diabatization = loadDiabatizationSettings(jss["diabatization"])
  utilitySettings = loadUtilitySettings(jss["utility"])
  interpolationSettings = loadInterpolationSettings(jss["interpolation"])

  settings = CalculationSettings(strategy_name,
    coordinate_start, coordinate_step, coordinate_piece, coordinate_step_error,
    asymptoticSettings, nonadiabaticAreas, diabatization, utilitySettings, interpolationSettings)
  return settings
end

function loadInterpolationSettings(js)
  interpolation_H = get(MAPPING_InterpolationType, js["hamiltonian"], UNSUPPORTED_INTERP::InterpolationType)
  interpolation_∂_∂R = get(MAPPING_InterpolationType, js["coupling-∂_∂R"], UNSUPPORTED_INTERP::InterpolationType)
  if UNSUPPORTED_INTERP::InterpolationType == interpolation_H
    throw(DomainError("Unsupported hamiltonian interpolation sort: $(js["hamiltonian"])"))
  end
  if UNSUPPORTED_INTERP::InterpolationType == interpolation_∂_∂R
    throw(DomainError("Unsupported ∂_∂R-coupling interpolation sort: $(js["coupling-∂_∂R"])"))
  end
  return InterpolationSettings(interpolation_H, interpolation_∂_∂R)
end

function loadUtilitySettings(js)
  channel_ionic_number = js["channel-ionic-number"]
  channel_lowest_number = js["channel-lowest-number"]
  return UtilitySettings(channel_ionic_number, channel_lowest_number)
end

function loadDiabatizationSettings(js)
  step_min = abs(js["coordinate-step"]["min"])
  step_max = abs(js["coordinate-step"]["max"])
  @assert step_min <= step_max
  use_prev_expression = haskey(js, "use-last-transformation-matrix-from") ?
    Nullable{Float64}(js["use-last-transformation-matrix-from"]) : Nullable{Float64}()

  js_selected_areas = js["selected-areas"]
  selected_areas = Vector{SelectedDiabatizationArea}(undef, 0)
  for area in js_selected_areas
      state1 = area["states"][1]
      state2 = area["states"][2]
      extra_length = area["extra-length"]
      diab_area = SelectedDiabatizationArea(
        area["coordinate"], (state1, state2), extra_length, area["bunch-exclude"])
      push!(selected_areas, diab_area)
  end
  sort!(selected_areas, by = diab_area -> diab_area.coordinate)

  return DiabatizationSettings(
    selected_areas,
    js["keep-initial-conditions"],
    js["coordinate-start"],
    (step_min, step_max),
    js["coordinate-stop"],
    js["area-closeness"],
    use_prev_expression
  )
end

function loadNonadiabaticAreas(js, p_coordinate_start, p_coordinate_step, p_coordinate_piece, p_coordinate_step_error)
  if (!haskey(js, "coordinate-start") && p_coordinate_start == -1) || (js["coordinate-start"] < 0 && p_coordinate_start == -1)
    throw(ErrorException("Please define the setting 'coordinate-start' in the global 'settings' section or in the 'nonadiabatic-areas' child configuration element using a small positive number."))
  end
  if (!haskey(js, "coordinate-step") && p_coordinate_step == -1)  || (js["coordinate-step"] < 0 && p_coordinate_step == -1)
    throw(ErrorException("Please define the setting 'coordinate-step' in the global 'settings' section or in the 'nonadiabatic-areas' child configuration element using a small positive number."))
  end
  if (!haskey(js, "coordinate-piece") && p_coordinate_step == -1)  || (js["coordinate-piece"] < 0 && p_coordinate_step == -1)
    throw(ErrorException("Please define the setting 'coordinate-piece' in the global 'settings' section or in the 'nonadiabatic-areas' child configuration element using a small positive number."))
  end
  if (!haskey(js, "coordinate-step-error") && p_coordinate_step_error == -1)  || (js["coordinate-step-error"] < 0 && p_coordinate_step_error == -1)
    throw(ErrorException("Please define the setting 'coordinate-step-error' in the global 'settings' section or in the 'nonadiabatic-areas' child configuration element using a small positive number."))
  end

  coordinate_start = deriveCoordinateParameter(js, "coordinate-start", p_coordinate_start)
  coordinate_step = deriveCoordinateParameter(js, "coordinate-step", p_coordinate_step)
  coordinate_piece = deriveCoordinateParameter(js, "coordinate-piece", p_coordinate_piece)
  coordinate_step_error = deriveCoordinateParameter(js, "coordinate-step-error", p_coordinate_step_error)
  if coordinate_step < 1e-9
    throw(ErrorException("The configured step by coordinate ΔR=$coordinate_step is too small or negative."))
  end
  if coordinate_piece <= coordinate_step
    throw(ErrorException("The configured coordidate piece $coordinate_piece is less than ΔR=$coordinate_step."))
  end
  if coordinate_step_error < 1e-12
    throw(ErrorException("The configured step error by coordinate ΔR=$coordinate_step_error is too small or negative."))
  end

  js_areas = js["areas"]
  areas = Dict{NonadiabaticAreaTypes, NonadiabaticAreaSettings}()
  for area in js_areas
    areaType = get(MAPPING_NonadiabaticAreaTypes, area.first, UNSUPPORTED::NonadiabaticAreaTypes)
    if areaType == UNSUPPORTED::NonadiabaticAreaTypes
      throw(ErrorException("Unsupported non-adiabatic area type '$(area.first)'"))
    end
    if areaType == SINGLE_PEAK::NonadiabaticAreaTypes
      areaSettings = buildSinglePeakNonadiabaticArea(area.second)
      areas[areaType] = areaSettings
    end
  end

  areasConfig = NonadiabaticAreasConfiguration(coordinate_start, coordinate_step, coordinate_piece, coordinate_step_error, areas)
  return areasConfig
end

function buildSinglePeakNonadiabaticArea(area)
  error_∂_∂R_peak = area["error-∂_∂R-peak"]
  vanishing_∂_∂R_value = area["vanishing-∂_∂R-value"]
  error_vanishing_∂_∂R_value = area["error-vanishing-∂_∂R-value"]
  error_potential_distance_minimal = area["error-potential-distance-minimal"]
  error_potential_distance_coordinate = area["error-potential-distance-coordinate"]
  error_potential_∂_∂R_coordinate = area["error-potential-∂_∂R-coordinate"]
  return SinglePeakNonadiabaticAreaSettings(
    error_∂_∂R_peak,
    vanishing_∂_∂R_value,
    error_vanishing_∂_∂R_value,
    error_potential_distance_minimal,
    error_potential_distance_coordinate,
    error_potential_∂_∂R_coordinate
  )
end

function loadAsymptotics(js, p_coordinate_start, p_coordinate_step, p_coordinate_step_error)
  if (!haskey(js, "coordinate-start") && p_coordinate_start == -1) || (js["coordinate-start"] < 0 && p_coordinate_start == -1)
    throw(DomainError("Please define the setting 'coordinate-start' in the global 'settings' section or in the 'asymptotics' child configuration element using a small positive number."))
  end
  if (!haskey(js, "coordinate-step") && p_coordinate_step == -1)  || (js["coordinate-step"] < 0 && p_coordinate_start == -1)
    throw(DomainError("Please define the setting 'coordinate-step' in the global 'settings' section or in the 'asymptotics' child configuration element using a small positive number."))
  end
  if (!haskey(js, "coordinate-step-error") && p_coordinate_step_error == -1)  || (js["coordinate-step-error"] < 0 && p_coordinate_step_error == -1)
    throw(DomainError("Please define the setting 'coordinate-step-error' in the global 'settings' section or in the 'nonadiabatic-areas' child configuration element using a small positive number."))
  end

  coordinate_start = deriveCoordinateParameter(js, "coordinate-start", p_coordinate_start)
  coordinate_step = deriveCoordinateParameter(js, "coordinate-step", p_coordinate_step)
  coordinate_step_error = deriveCoordinateParameter(js, "coordinate-step-error", p_coordinate_step_error)
  if coordinate_step < 1e-8
    throw(DomainError("The configured step by coordinate ΔR=$coordinate_step is too small or negative."))
  end
  if coordinate_step_error < 1e-12
    throw(DomainError("The configured step error by coordinate ΔR=$coordinate_step_error is too small or negative."))
  end

  coordinate_safety_step = js["coordinate-safety-step"]
  potential_asymptotic_value_error = js["potential-asymptotic-value-error"]
  ∂_∂R_asymptotic_value_error = js["∂_∂R-asymptotic-value-error"]
  ∂_∂R_zero_value_error = js["∂_∂R-zero-value-error"]
  return AsymptoticSettings(
    coordinate_start,
    coordinate_step,
    coordinate_step_error,
    coordinate_safety_step,
    potential_asymptotic_value_error,
    ∂_∂R_asymptotic_value_error,
    ∂_∂R_zero_value_error
  )
end

function deriveCoordinateParameter(js, name::AbstractString, parent_value)
  if haskey(js, name)
    value = js[name]
    if value > 0 && parent_value < 0
      return value
    elseif value < 0 && parent_value > 0
      return parent_value
    elseif value > 0 && parent_value > 0
      return value
    else
      throw(ErrorException("Unexpected behaviour; name=$name, parent_value=$parent_value."))
    end
  elseif !haskey(js, name) && parent_value > 0
    return parent_value
  else
    throw(ErrorException("Unexpected behaviour; name=$name, parent_value=$parent_value."))
  end
end

function deriveStrategy(strategyName::AbstractString)
  strategy = get(MAPPING_CalculationStrategy, strategyName, UNKNOWN::CalculationStrategy)
  if strategy == UNKNOWN::CalculationStrategy
    throw(DomainError("Illegal strategy name '$strategyName'."))
  end
  return strategy
end

function loadRawData(filePath::AbstractString)
  return DataFrame(load(filePath;
    delim=' ', spacedelim=true, header_exists=true))
end
