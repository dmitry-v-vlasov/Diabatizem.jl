# ---------------------------------
# Miscellaneous utility functions
# ---------------------------------
"""
  Calculates a column number in a data file for with a given
  symmetric or antisymmetric (square) matrix with ZEROs in the main diagonal by
  its line an column indices.

  * i - matrix line
  * j - matrix column
  * N - matrix size

  returns: see sescription
  throws:
    * DomainError in case i is equal to j
"""
function dataColumnOfSymetricMatrix(i::Int, j::Int, N::Int)
  if i < j
    return (2*(N - 1) - i) * (i - 1) / 2 + (j - 1)
  elseif i > j
    return (2*(N - 1) - j) * (j - 1) / 2 + (i - 1)
  else
    throw(DomainError("The equal line and column value '$i' is undefined for the (anti)symmetric matrix with zero main diagonal."))
  end
end

"""
Convert vector element number to a pair of N×N matrix indices.
"""
function mpos(l::Int, N::Int)
  dv = divrem(l, N)
  return dv[1] == 0 ? (1, dv[2]) : (dv[1] + 1, dv[2])
end

"""
Convert a pair of N×N matrix indices to a one-dimentional array index.
"""
function mvec(i::Int, j::Int, N::Int)
  return N*(i - 1) + j
end
