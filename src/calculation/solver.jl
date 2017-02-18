using Sundials

function diabatize(Hₐ::Array{Function, 2}, ∂_∂R::Array{Function, 2}, ∂_∂Rᵐᵒᵈᵉˡ::Array{Array{LandauZenerArea, 1}, 2}, config::DiabatizationSettings)
  Rᵇᵉᵍⁱⁿ = config.coordinate_start; Rᵉⁿᵈ = config.coordinate_stop
  ΔRᵐᵃˣ = config.coordinate_step

end
