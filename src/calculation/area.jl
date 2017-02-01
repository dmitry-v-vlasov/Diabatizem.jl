using Calculus

function detectSinglePeakAreas(M_∂_∂R::Array{Function, 2}, settings::SinglePeakNonadiabaticAreaSettings, Rstop::Float64)
  ϵ_peak = settings.error_∂_∂R_peak
  N = size(M_∂_∂R, 1)
  for i = 1:N, j = 1:N
    if i < j && j - i == 1
      points = Float64[]
      τ = M_∂_∂R[i, j]
      Rₛ = 0.1; ΔR = 0.001
      τᵥ = τ(Rₛ); ∂τ_∂R = derivative(τ, Rₛ); ∂²τ_∂R² = second_derivative(τ, Rₛ)
      for R = Rₛ+ΔR:ΔR:Rstop
        τᵥ_ΔR = τ(R); ∂τ_∂R_ΔR = derivative(τ, R); ∂²τ_∂R²_ΔR = second_derivative(τ, R)
        # -----------

        sgn_∂τ_∂R = sign(∂τ_∂R); sgn_∂τ_∂R_ΔR = sign(∂τ_∂R_ΔR)
        #sgn_∂²τ_∂R² = sign(∂²τ_∂R²); sgn_∂²τ_∂R²_ΔR = sign(∂²τ_∂R²_ΔR)

        Δτ = abs(τᵥ_ΔR - τᵥ)
        if sgn_∂τ_∂R == -sgn_∂τ_∂R_ΔR
          if Δτ < ϵ_peak
            push!(points, (R + (R + ΔR))/2)
            R += ΔR
          else
            R -= ΔR
            ΔR /= 2
          end
        else
          if Δτ < ϵ_peak
            ΔR *= 2
          end
        end

        # -----------
        τᵥ = τᵥ_ΔR; ∂τ_∂R = ∂τ_∂R_ΔR; ∂²τ_∂R² = ∂²τ_∂R²_ΔR
      end
    end
  end
end
