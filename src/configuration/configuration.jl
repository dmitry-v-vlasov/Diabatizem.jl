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

type UtilitySettings
  channel_ionic_number::Int
  channel_lowest_number::Int
end

type AsymptoticSettings
  coordinate_start::Float64
  coordinate_safety_step::Float64
  potential_asymptotic_value_error::Float64
  ∂_∂R_asymptotic_value_error::Float64
  ∂_∂R_zero_value_error::Float64
end

abstract NonadiabaticAreaSettings
type SinglePeakNonadiabaticAreaSettings <: NonadiabaticAreaSettings
  error_∂_∂R_peak::Float64
  vanishing_∂_∂R_value::Float64
  error_potential_distance_minimal::Float64
  error_potential_distance_coordinate::Float64
  error_potential__∂_∂R_coordinate::Float64
end

type CalculationSettings
  strategy::CalculationStrategy
  asymptotics::AsymptoticSettings
  nonadiabatic_areas::NonadiabaticAreaSettings
  utility::UtilitySettings
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
  input_data = InputData(loadRawData(input_paths.file_hamiltonian_adiabatic), loadRawData(input_paths.file_coupling_∂_∂R_adiabatic), Nullable{DataFrame}())

  output_paths = OutputPaths(js["output-data"]["hamiltonian-diabatic"], js["output-data"]["coupling-∂_∂R-diabatic"], js["output-data"]["transformation-matrix"])

  strategy_name = CalculationStrategy(deriveStrategy(js["settings"]["strategy"]))
  if LANDAU_ZENER_WITH_EXTERNAL_MODEL_DATA == strategy_name || EXTERNAL_MODEL_DATA == strategy_name
    if haskey(js["input-data"], "coupling_∂_∂R_adiabatic_model")
      input_paths.file_coupling_∂_∂R_adiabatic_model = js["input-data"]["coupling_∂_∂R_adiabatic_model"]
    else
      throw(DomainError("The configuration setting 'coupling_∂_∂R_adiabatic_model' is obligatory for the strategy $strategy_name"))
    end
  end
  if LANDAU_ZENER_WITH_EXTERNAL_MODEL_DATA == strategy_name || EXTERNAL_MODEL_DATA == strategy_names
    input_data.coupling_∂_∂R_adiabatic_model = loadRawData(get(input_paths.file_coupling_∂_∂R_adiabatic_model))
  end

  settings = CalculationSettings(strategy_name)

  return Configuration(input_paths, input_data, output_paths, settings)
end

function deriveStrategy(strategyName::AbstractString)
  strategy = begin
    if "Landau-Zener" == strategyName
      return LANDAU_ZENER
    elseif "Landau-Zener-with-external-model-data" == strategyName
      return LANDAU_ZENER_WITH_EXTERNAL_MODEL_DATA
    elseif "external-model-data" == strategyName
      return EXTERNAL_MODEL_DATA
    elseif "just-make-them-zeros" == strategyName
      return JUST_MAKE_THEM_ZEROS
    else
      throw(DomainError("Illegal strategy name '$strategyName'."))
    end
  end
end

function loadRawData(filePath::AbstractString)
  return DataFrames.readtable(filePath, header = true, separator = ' ', allowcomments = true, commentmark = '#', skipblanks = true, encoding = :utf8, normalizenames = true)
end
