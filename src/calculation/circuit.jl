# -----------------------------------------------------------------------------
# Derivation of the initial parameters for the calculation schema
# -----------------------------------------------------------------------------

using Calculus
using Formatting

function potentialAsymptoticValues(potentials::Array{Function, 2},
                                    coordinate_start::Float64,
                                    coordinate_step::Float64,
                                    coordinate_stop::Float64,
                                    safety_coordinate_step::Float64,
                                    asymptotic_value_error::Float64)
  N = size(potentials, 1)
  potentialAsymptotics = Vector{PotentialAsymptotic}(undef, N)

  R₀ = coordinate_start
  ΔR = coordinate_step
  Rₐ = coordinate_stop
  Δh = safety_coordinate_step
  ϵ = abs(asymptotic_value_error)
  for l = 1:N
    U = potentials[l][l]
    zero_counter = 0
    U_R_prev = U(R₀-ΔR)
    for R = R₀:ΔR:Rₐ
      U_R = U(R)
      ΔU = U_R - U_R_prev
      dU_dR = derivative(U, R); d²U_dR² = second_derivative(U, R)

      if abs(ΔU) < ϵ && abs(dU_dR) < ϵ && abs(d²U_dR²) < ϵ
        zero_counter += 1
        if zero_counter == 10
          ΔU_Δh = U(R + Δh) - U_R
          dU_dR_Δh = derivative(U, R + Δh); d²U_dR²_Δh = second_derivative(U, R + Δh)
          if abs(ΔU_Δh) < ϵ && abs(dU_dR_Δh) < ϵ && abs(d²U_dR²_Δh) < ϵ
            potentialAsymptotics[l] = PotentialAsymptotic(l, U_R, R)
            break
          end
        end
      end

      if R == Rₐ && zero_counter < 10
        error("The potential with the number $l (possible state name is '$(⚛⚛_STATES[l])') has unclear asymptotic behaviour.")
      end
    end
  end

  return potentialAsymptotics
end

function toStringPotentialAsymptotics(asymptotics::Array{PotentialAsymptotic})
  buff = IOBuffer()
  for asymptotic in asymptotics
    ψ = asymptotic.state
    U_∞ = sprintf1("%.4f", asymptotic.value)
    Rₐ = sprintf1("%.4f", asymptotic.coordinate)
    print(buff, "{$(ψ)→$(⚛⚛_STATES[ψ]): U(∞)=$(U_∞); Rₐ=$(Rₐ)}")
    print(buff, "\n")
  end
  return takebuf_string(buff)
end
