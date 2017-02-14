using JSON
using DataFrames
using DataStructures

const EMPTY = ""

type InputPaths
  file_hamiltonian_adiabatic::AbstractString
  file_coupling_∂_∂R_adiabatic::AbstractString
  file_coupling_∂_∂R_adiabatic_model::Nullable{AbstractString}
end

type InputData
  hamiltonian_adiabatic::DataFrame
  coupling_∂_∂R_adiabatic::DataFrame
  coupling_∂_∂R_adiabatic_model::Nullable{DataFrame}
end

type OutputPaths
  file_hamiltonian_diabatic::AbstractString
  file_coupling_∂_∂R_diabatic::AbstractString
  file_transformation_matrix::AbstractString
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

type UtilitySettings
  channel_ionic_number::Int
  channel_lowest_number::Int
end

@enum InterpolationType SPLINE_LINEAR SPLINE_QUADRATIC SPLINE_CUBIC UNSUPPORTED_INTERP
const MAPPING_InterpolationType = Dict(
  "spline-linear" => SPLINE_LINEAR::InterpolationType,
  "spline-quadratic" => SPLINE_QUADRATIC::InterpolationType,
  "spline-cubic" => SPLINE_CUBIC::InterpolationType
)
type InterpolationSettings
  hamiltonian::InterpolationType
  coupling_∂_∂R::InterpolationType
end

type AsymptoticSettings
  coordinate_start::Float64
  coordinate_step::Float64
  coordinate_step_error::Float64
  coordinate_safety_step::Float64
  potential_asymptotic_value_error::Float64
  ∂_∂R_asymptotic_value_error::Float64
  ∂_∂R_zero_value_error::Float64
end

abstract NonadiabaticAreaSettings
type SinglePeakNonadiabaticAreaSettings <: NonadiabaticAreaSettings
  error_∂_∂R_peak::Float64
  vanishing_∂_∂R_value::Float64
  error_vanishing_∂_∂R_value::Float64
  error_potential_distance_minimal::Float64
  error_potential_distance_coordinate::Float64
  error_potential_∂_∂R_coordinate::Float64
end
type NonadiabaticAreasConfiguration
  coordinate_start::Float64
  coordinate_step::Float64
  coordinate_step_error::Float64
  nonadiabatic_areas::Dict{NonadiabaticAreaTypes, NonadiabaticAreaSettings}
end

type CalculationSettings
  strategy::CalculationStrategy
  coordinate_start::Float64
  coordinate_step::Float64
  coordinate_step_error::Float64
  asymptotics::AsymptoticSettings
  nonadiabatic_areas::NonadiabaticAreasConfiguration
  utility::UtilitySettings
  interpolation::InterpolationSettings
end

type Configuration
  input_paths::InputPaths
  input_data::InputData
  output_paths::OutputPaths
  settings::CalculationSettings
end

function loadConfiguration(filePath::AbstractString)
  js = JSON.parsefile(filePath; dicttype=DataStructures.OrderedDict)

  input_paths = InputPaths(js["input-data"]["hamiltonian-adiabatic"], js["input-data"]["coupling_∂_∂R_adiabatic"], Nullable{AbstractString}())

  Hₐ_data = loadRawData(input_paths.file_hamiltonian_adiabatic)
  fixPotentialAsymptotics!(Hₐ_data, loadUtilitySettings(js["settings"]["utility"]))
  input_data = InputData(Hₐ_data, loadRawData(input_paths.file_coupling_∂_∂R_adiabatic), Nullable{DataFrame}())

  output_paths = OutputPaths(js["output-data"]["hamiltonian-diabatic"], js["output-data"]["coupling-∂_∂R-diabatic"], js["output-data"]["transformation-matrix"])

  settings = loadCalculationSettings(js, input_paths, input_data)

  return Configuration(input_paths, input_data, output_paths, settings)
end

function fixPotentialAsymptotics!(Hₐ_data::DataFrame, utilitySettings)
  Nₚ = size(Hₐ_data, 1)
  Nᵩ = size(Hₐ_data, 2) - 1
  lowest_channel = utilitySettings.channel_lowest_number > 0 ? utilitySettings.channel_lowest_number : 1
  Uₗ_∞ = Hₐ_data[Nₚ, lowest_channel + 1]
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
  if LANDAU_ZENER_WITH_EXTERNAL_MODEL_DATA::CalculationStrategy == strategy_name || EXTERNAL_MODEL_DATA::CalculationStrategy == strategy_names
    input_data.coupling_∂_∂R_adiabatic_model = loadRawData(get(input_paths.file_coupling_∂_∂R_adiabatic_model))
  end

  coordinate_start = haskey(jss, "coordinate-start") ? (jss["coordinate-start"] < 0 ? -1 : jss["coordinate-start"]) : -1
  coordinate_step = haskey(jss, "coordinate-step") ? (jss["coordinate-step"] < 0 ? -1 : jss["coordinate-step"]) : -1
  coordinate_step_error = haskey(jss, "coordinate-step-error") ? (jss["coordinate-step-error"] < 0 ? -1 : jss["coordinate-step-error"]) : -1

  asymptoticSettings = loadAsymptotics(jss["asymptotics"], coordinate_start, coordinate_step, coordinate_step_error)
  nonadiabaticAreas = loadNonadiabaticAreas(jss["nonadiabatic-areas"], coordinate_start, coordinate_step, coordinate_step_error)
  utilitySettings = loadUtilitySettings(jss["utility"])
  interpolationSettings = loadInterpolationSettings(jss["interpolation"])

  settings = CalculationSettings(strategy_name,
    coordinate_start, coordinate_step, coordinate_step_error,
    asymptoticSettings, nonadiabaticAreas, utilitySettings, interpolationSettings)
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

function loadNonadiabaticAreas(js, p_coordinate_start, p_coordinate_step, p_coordinate_step_error)
  if (!haskey(js, "coordinate-start") && p_coordinate_start == -1) || (js["coordinate-start"] < 0 && p_coordinate_start == -1)
    throw(DomainError("Please define the setting 'coordinate-start' in the global 'settings' section or in the 'nonadiabatic-areas' child configuration element using a small positive number."))
  end
  if (!haskey(js, "coordinate-step") && p_coordinate_step == -1)  || (js["coordinate-step"] < 0 && p_coordinate_step == -1)
    throw(DomainError("Please define the setting 'coordinate-step' in the global 'settings' section or in the 'nonadiabatic-areas' child configuration element using a small positive number."))
  end
  if (!haskey(js, "coordinate-step-error") && p_coordinate_step_error == -1)  || (js["coordinate-step-error"] < 0 && p_coordinate_step_error == -1)
    throw(DomainError("Please define the setting 'coordinate-step-error' in the global 'settings' section or in the 'nonadiabatic-areas' child configuration element using a small positive number."))
  end

  coordinate_start = deriveCoordinateParameter(js, "coordinate-start", p_coordinate_start)
  coordinate_step = deriveCoordinateParameter(js, "coordinate-step", p_coordinate_step)
  coordinate_step_error = deriveCoordinateParameter(js, "coordinate-step-error", p_coordinate_step_error)
  if coordinate_step < 1e-9
    throw(DomainError("The configured step by coordinate ΔR=$coordinate_step is too small or negative."))
  end
  if coordinate_step_error < 1e-12
    throw(DomainError("The configured step error by coordinate ΔR=$coordinate_step_error is too small or negative."))
  end

  js_areas = js["areas"]
  areas = Dict{NonadiabaticAreaTypes, NonadiabaticAreaSettings}()
  for area in js_areas
    areaType = get(MAPPING_NonadiabaticAreaTypes, area.first, UNSUPPORTED::NonadiabaticAreaTypes)
    if areaType == UNSUPPORTED::NonadiabaticAreaTypes
      throw(DomainError("Unsupported non-adiabatic area type '$(area.first)'"))
    end
    if areaType == SINGLE_PEAK::NonadiabaticAreaTypes
      areaSettings = buildSinglePeakNonadiabaticArea(area.second)
      areas[areaType] = areaSettings
    end
  end

  areasConfig = NonadiabaticAreasConfiguration(coordinate_start, coordinate_step, coordinate_step_error, areas)
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
  return DataFrames.readtable(
    filePath,
    header = true, separator = ' ',
    allowcomments = true, commentmark = '#',
    skipblanks = true, encoding = :utf8, normalizenames = true
  )
end
