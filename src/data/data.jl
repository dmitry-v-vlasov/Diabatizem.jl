using DataFrames
using Calculus
using Formatting

import Dierckx

type Data
  Hₐ::Array{Function, 2}
  ∂_∂R::Array{Function, 2}
  itp_Hₐ::Array{Dierckx.Spline1D, 2}
  itp_∂_∂R::Array{Dierckx.Spline1D, 2}
  function Data() new() end
end

function saveData(Rᵖᵒⁱⁿᵗˢ::Vector{Float64},
    S::Vector{Array{Float64, 2}},
    Sᵈᵃᵗᵃ::Array{Float64, 2},
    Uᴰᵈᵃᵗᵃ::Array{Float64, 2},
    Hᴰᵈᵃᵗᵃ::Array{Float64, 2},
    ∂_∂Rᴰᵈᵃᵗᵃ::Array{Float64, 2},
    ∂²_∂R²ᴰᵈᵃᵗᵃ::Array{Float64, 2},
    ∂²_∂R²ᴰᵈᵃᵗᵃ_diag::Array{Float64, 2},
    out::OutputPaths)
  Logging.configure(level=INFO)

  @assert length(Rᵖᵒⁱⁿᵗˢ) == size(Uᴰᵈᵃᵗᵃ, 1)
  @assert length(Rᵖᵒⁱⁿᵗˢ) == size(Hᴰᵈᵃᵗᵃ, 1)
  @assert length(Rᵖᵒⁱⁿᵗˢ) == size(∂_∂Rᴰᵈᵃᵗᵃ, 1)
  @assert length(Rᵖᵒⁱⁿᵗˢ) == size(∂²_∂R²ᴰᵈᵃᵗᵃ, 1)
  @assert length(Rᵖᵒⁱⁿᵗˢ) == size(∂²_∂R²ᴰᵈᵃᵗᵃ_diag, 1)

  @assert size(Hᴰᵈᵃᵗᵃ, 2) == size(Uᴰᵈᵃᵗᵃ, 2)*(size(Uᴰᵈᵃᵗᵃ, 2) - 1)/2
  @assert size(Hᴰᵈᵃᵗᵃ, 2) == size(∂_∂Rᴰᵈᵃᵗᵃ, 2)
  @assert size(Hᴰᵈᵃᵗᵃ, 2) == size(∂²_∂R²ᴰᵈᵃᵗᵃ, 2)
  @assert size(Uᴰᵈᵃᵗᵃ, 2) == size(∂²_∂R²ᴰᵈᵃᵗᵃ_diag, 2)

  L = length(Rᵖᵒⁱⁿᵗˢ)
  N = size(Uᴰᵈᵃᵗᵃ, 2)

  # -----------
  Sᵗᵃᵇˡᵉ = makeMatrixElementTable(Rᵖᵒⁱⁿᵗˢ, Sᵈᵃᵗᵃ, :general, "C", N)
  # -----------
  Hᵈⁱᵃᵍ = makeMatrixElementTable(Rᵖᵒⁱⁿᵗˢ, Uᴰᵈᵃᵗᵃ, :diagonal, "H", N)
  # -----------
  Hᵒᶠᶠᵈⁱᵃᵍ = makeMatrixElementTable(Rᵖᵒⁱⁿᵗˢ, Hᴰᵈᵃᵗᵃ, :symmetric, "H", N)
  # -----------
  ∂_∂R = makeMatrixElementTable(Rᵖᵒⁱⁿᵗˢ, ∂_∂Rᴰᵈᵃᵗᵃ, :antisymmetric, "d/dR", N)
  # -----------
  ∂²_∂R² = makeMatrixElementTable(Rᵖᵒⁱⁿᵗˢ, ∂²_∂R²ᴰᵈᵃᵗᵃ, :symmetric, "d2/dR2", N)
  # -----------
  ∂²_∂R²ᵈⁱᵃᵍ = makeMatrixElementTable(Rᵖᵒⁱⁿᵗˢ, ∂²_∂R²ᴰᵈᵃᵗᵃ_diag, :diagonal, "d2/dR2", N)
  # -----------

  # -----------
  info("Saving the transformation matrix to '$(out.file_transformation_matrix_general)'")
  saveMatrixList(Rᵖᵒⁱⁿᵗˢ, S, out.file_transformation_matrix_general)
  info("Saving the transformation matrix table to '$(out.file_transformation_matrix)'")
  saveMatrixElementTable(Sᵗᵃᵇˡᵉ, out.file_transformation_matrix)
  info("Saving the transformation matrix table to '$(out.file_potentials_diabatic)'")
  saveMatrixElementTable(Hᵈⁱᵃᵍ, out.file_potentials_diabatic)
  info("Saving the transformation matrix table to '$(out.file_hamiltonian_diabatic)'")
  saveMatrixElementTable(Hᵒᶠᶠᵈⁱᵃᵍ, out.file_hamiltonian_diabatic)
  info("Saving the transformation matrix table to '$(out.file_coupling_∂_∂R_diabatic)'")
  saveMatrixElementTable(∂_∂R, out.file_coupling_∂_∂R_diabatic)
  info("Saving the transformation matrix table to '$(out.file_coupling_∂²_∂R²_diabatic)'")
  saveMatrixElementTable(∂²_∂R², out.file_coupling_∂²_∂R²_diabatic)
  info("Saving the transformation matrix table to '$(out.file_coupling_∂²_∂R²_diabatic_diag)'")
  saveMatrixElementTable(∂²_∂R²ᵈⁱᵃᵍ, out.file_coupling_∂²_∂R²_diabatic_diag)
  info("Done")
  # -----------
end

function saveMatrixList(Rᵖᵒⁱⁿᵗˢ::Vector{Float64}, M::Vector{Array{Float64, 2}}, file_name::AbstractString)
  @assert length(Rᵖᵒⁱⁿᵗˢ) == length(M) "$(length(Rᵖᵒⁱⁿᵗˢ))≠$(length(M))"
  fos = open(file_name, "a")
  L = length(Rᵖᵒⁱⁿᵗˢ)
  for l = 1:L
    Rˡ = Rᵖᵒⁱⁿᵗˢ[l]; Mˡ = M[l]
    N = size(Mˡ, 1)
    for i = 1:N, j = 1:N
      if abs(Mˡ[i, j]) < 1e-15
        Mˡ[i, j] = 0.0
      elseif abs(1-abs(Mˡ[i, j])) < 1e-5
        Mˡ[i, j] = sign(Mˡ[i, j]) * 1.0
      else
        Mˡ[i, j] = Mˡ[i, j]
      end
    end

    fmtrfunc = generate_formatter("%11.9e")
    buf = IOBuffer()
    for i = 1:N, j = 1:N
      print(buf, fmtrfunc(Mˡ[i, j]))
      if i <= N && j < N
        print(buf, "  ")
      elseif i <= N && j == N
        print(buf, "\n")
      end
    end

    write(fos, "=======================================================\n")
    write(fos, "R = $Rˡ\n")
    write(fos, "-----------\n")
    write(fos, takebuf_string(buf))
    write(fos, "\n")
  end
  flush(fos)
  close(fos)
end

function saveMatrixElementTable(data::DataFrame, file_name::AbstractString)
    writetable(file_name, data; separator=' ', quotemark=' ', header=true, nastring="EMPTY")
end

function makeMatrixElementTable(Rᵖᵒⁱⁿᵗˢ::Vector{Float64}, A::Array{Float64, 2}, dataType::Symbol, operator::AbstractString, N::Int)
  @assert length(Rᵖᵒⁱⁿᵗˢ) == size(A, 1)
  data = DataFrame()
  data[:R] = Rᵖᵒⁱⁿᵗˢ
  Nᶜᵒˡ = size(A, 2)
  if dataType == :diagonal
    @assert Nᶜᵒˡ == N "$Nᶜᵒˡ≠$N"
    for l = 1:Nᶜᵒˡ
      state = ⚛⚛_STATES[l]
      data[Symbol("<$state|$operator|$state>")] = A[:, l]
    end
  elseif (dataType == :symmetric || dataType == :antisymmetric)
    @assert Nᶜᵒˡ == N*(N-1)/2 "$(Nᶜᵒˡ)≠$(N*(N-1)/2), N=$N"
    lᵖ = 1
    for i = 1:N, j=i+1:N
      bra = ⚛⚛_STATES[i]; ket = ⚛⚛_STATES[j]
      l = dataColumnOfSymetricMatrix(i, j, N)
      @assert lᵖ <= l "$lᵖ>$l"; lᵖ = l;
      @assert l <= Nᶜᵒˡ "$l>$Nᶜᵒˡ"
      data[Symbol("<$bra|$operator|$ket>")] = A[:, l]
    end
  elseif dataType == :general
    @assert Nᶜᵒˡ == N*N "$(Nᶜᵒˡ)≠$(N*N), N=$N"
    for l = 1:Nᶜᵒˡ
      i, j = mpos(l, N)
      bra = ⚛⚛_STATES[i]; ket = ⚛⚛_STATES[j]
      data[Symbol("<$bra|$operator|$ket>")] = A[:, l]
    end
  else
    error("Unsupported data type.")
  end
  return data
end

function buildData(table_Hₐ::DataFrame, table_∂_∂R::DataFrame, interpolationSettings::InterpolationSettings)
  data = Data()
  buildHₐ!(table_Hₐ, data, interpolationSettings.hamiltonian)
  build_d_dR!(table_∂_∂R, data, numberOfChannels(table_Hₐ), interpolationSettings.coupling_∂_∂R)
  return data
end

function buildHₐ!(table_Hₐ::DataFrame, data::Data, interpolationType::InterpolationType)
  Logging.configure(level=INFO)
  info("Doing spline interpolation for adiabatic hamiltonian")
  X = convert(Vector{Float64}, table_Hₐ[1])
  ΔR = X[2] - X[1]
  N = numberOfChannels(table_Hₐ)
  data.Hₐ = Array{Function, 2}(N, N)
  data.itp_Hₐ = Array{Dierckx.Spline1D, 2}(N, N)
  for i = 1:N, j = 1:N
    if i == j
      Y = convert(Vector{Float64}, table_Hₐ[i + 1])
      spl = Dierckx.Spline1D(X, Y; w=ones(length(X)), k=splineDegree(interpolationType), bc="nearest", s=0.0)
      setindex!(data.Hₐ, R -> Dierckx.evaluate(spl, R), i, j)
      setindex!(data.itp_Hₐ, spl, i, j)
      info("Spline initialized; i=$i, j=$j. Asymptotics: table→ ($(X[end]), $(Y[end])), spline→ ($(X[end]), $(Dierckx.evaluate(spl, X[end])))")
    else
      setindex!(data.Hₐ, R -> 0, i, j)
    end
  end
end

function build_d_dR!(table_∂_∂R::DataFrame, data::Data, N::Int, interpolationType::InterpolationType)
  Logging.configure(level=INFO)
  info("Doing spline interpolation for <|d/dR|> couplings")
  X = convert(Vector{Float64}, table_∂_∂R[1])
  ΔR = X[2] - X[1]
  Nc = size(X, 1) - 1
  data.∂_∂R = Array{Function, 2}(N, N)
  data.itp_∂_∂R = Array{Dierckx.Spline1D, 2}(N, N)
  for i = 1:N, j = 1:N
    if i ≠ j
      info("Interpolation for <$i|d/dR|$j>")
      l = dataColumnOfSymetricMatrix(i, j, N)
      Y = convert(Vector{Float64}, table_∂_∂R[l + 1])

      spl = Dierckx.Spline1D(X, (i < j ? Y : -Y); w=ones(length(X)), k=splineDegree(interpolationType), bc="nearest", s=0.0)
      setindex!(data.∂_∂R, R -> Dierckx.evaluate(spl, R), i, j)
      setindex!(data.itp_∂_∂R, spl, i, j)
    else
      setindex!(data.∂_∂R, R -> 0, i, j)
    end
  end
end

function loadInitialConditions(file_name::AbstractString, N::Int)
  if !isfile(file_name)
    return Nullable{Array{Float64, 2}}()
  end
  M = readdlm(file_name, ' ';
    header=false, skipstart=0, skipblanks=true,
    use_mmap=false, quotes=false, dims=(N, N),
    comments=true, comment_char='#')
  return Nullable{Array{Float64, 2}}(M)
end

function splineDegree(interpolationType::InterpolationType)
  if SPLINE_LINEAR::InterpolationType == interpolationType
    return 1
  elseif SPLINE_QUADRATIC::InterpolationType == interpolationType
    return 2
  elseif SPLINE_CUBIC::InterpolationType == interpolationType
    return 3
  end
end

function numberOfChannels(table_Hₐ::DataFrame)
  size(table_Hₐ, 2) - 1
end
