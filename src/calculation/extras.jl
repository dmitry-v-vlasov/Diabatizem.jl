import Dierckx

function calculate∂²_∂R²(Rᵖᵒⁱⁿᵗˢ::Vector{Float64},
    ∂_∂Rᴰᵈᵃᵗᵃ::Array{Float64, 2},
    N::Int)
    @assert length(Rᵖᵒⁱⁿᵗˢ) == size(∂_∂Rᴰᵈᵃᵗᵃ, 1)

    L = length(Rᵖᵒⁱⁿᵗˢ)
    M = size(∂_∂Rᴰᵈᵃᵗᵃ, 2)
    @assert M == dataSizeOfSymetricMatrix(N) "$M≠$(dataSizeOfSymetricMatrix(N))"

    ∂²_∂R²ᴰᵈᵃᵗᵃ = Array{Float64, 2}(L, M)
    ∂²_∂R²ᴰᵈᵃᵗᵃ_diag = Array{Float64, 2}(L, N)
    ∂_∂Rᴰᵈᵃᵗᵃ_func = Array{Function, 2}(N, N)
    ∂_∂Rᴰᵈᵃᵗᵃ_spl = Array{Dierckx.Spline1D, 2}(N, N)
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

    ∂²_∂R²ᴰᵈᵃᵗᵃ = Array{Float64, 2}(L, M)
    ∂²_∂R²ᴰᵈᵃᵗᵃ_diag = Array{Float64, 2}(L, N)
    for lʳ = 1:L
      τ⁽¹⁾ = matf2mat(Rᵖᵒⁱⁿᵗˢ[lʳ], ∂_∂Rᴰᵈᵃᵗᵃ_func)
      ∇τ⁽¹⁾ = matDerivative(Rᵖᵒⁱⁿᵗˢ[lʳ], ∂_∂Rᴰᵈᵃᵗᵃ_spl)
      τ⁽²⁾ = τ⁽¹⁾*τ⁽¹⁾ + ∇τ⁽¹⁾ # ???!!

      ∂²_∂R²ᴰᵈᵃᵗᵃ[lʳ, :] = matl2matupper(∇τ⁽¹⁾)
      ∂²_∂R²ᴰᵈᵃᵗᵃ_diag[lʳ, :] = diag(τ⁽²⁾)
    end

    return ∂²_∂R²ᴰᵈᵃᵗᵃ, ∂²_∂R²ᴰᵈᵃᵗᵃ_diag
end
