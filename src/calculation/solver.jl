using Sundials

function diabatize(Hₐ::Array{Function, 2}, ∂_∂R::Array{Function, 2}, ∂_∂Rᵐᵒᵈᵉˡ::Array{Array{LandauZenerArea, 1}, 2}, config::DiabatizationSettings)
  Rᵇᵉᵍⁱⁿ = config.coordinate_start; Rᵉⁿᵈ = config.coordinate_stop
  ΔRᵐᵃˣ = config.coordinate_step

end

function problemCauchy(;
  Y₀::Array{Float64} = nothing,
  Xᵖᵒⁱⁿᵗˢ::Vector{Float64} = nothing,
  prod_function::Function = nothing,
  ϵʳᵉˡ::Float64 = 1e-3,
  ϵᵃᵇˢ::Float64 = 1e-6)

  N = size(Y₀, 1)

  Yⁱⁿⁱᵗ = Vector{Float64}(N*N); fill!(Yⁱⁿⁱᵗ, 0)
  for i = 1:N, j = 1:N
    k = mvec(i, j); Yⁱⁿⁱᵗ[k] = Y₀[i, j]
  end

  Yʳᵉˢ = Sundials.cvode(prod_function, Yⁱⁿⁱᵗ, Xᵖᵒⁱⁿᵗˢ;
                          integrator = :BDF, reltol = ϵʳᵉˡ, abstol = ϵᵃᵇˢ)

  Nᵖᵒⁱⁿᵗˢ = size(Yʳᵉˢ, 1)
  Nˢᵒˡᵘᵗⁱᵒⁿˢ = size(Yʳᵉˢ, 2)
  if Nˢᵒˡᵘᵗⁱᵒⁿˢ ≠ N*N throw(ErrorException("Unexpected: Nˢᵒˡᵘᵗⁱᵒⁿˢ ≠ N×N: $Nˢᵒˡᵘᵗⁱᵒⁿˢ ≠ $N×$N ($(N*N))")) end
  Y = Array{Tuple{Float64, Array{Float64, Float64}}, 1}(Nᵖᵒⁱⁿᵗˢ)
  for k = 1:Nᵖᵒⁱⁿᵗˢ
    Yₖ = Array{Float64, Float64}(N*N)
    for l = 1:N*N
      (i, j) = mpos(l, N)
      Yₖ[i, j] = Yʳᵉˢ[l]
    end
    push!(Y, (Xᵖᵒⁱⁿᵗˢ[k], Yₖ))
  end

  return Y
end
