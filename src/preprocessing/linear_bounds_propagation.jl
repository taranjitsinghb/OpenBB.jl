# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-11T17:03:23+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: linear_bounds_propagation.jl
# @Last modified by:   massimo
# @Last modified time: 2019-07-05T17:16:11+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# single constraint linear bounds propagation
function bounds_propagation!(row::Int,
                            A::SparseMatrixCSC{Float64,Int},
                            cnsLoBs::Array{Float64,1},cnsUpBs::Array{Float64,1},
                            varLoBs::Array{Float64,1},varUpBs::Array{Float64,1},
                            dscIndices::Array{Int64,1})::Set{Int}
      cns = A[row,1:end]
      cnsLoB = cnsLoBs[row]
      cnsUpB = cnsUpBs[row]

      # use a set for the output
      updatedVars = Set{Int}()

      # get the sparsity of cns
      indices,coeffs = findnz(cns)

      # distiguish between positive and negative coefficients
      posCoeffs = Array{Float64,1}(undef,length(indices))
      @. posCoeffs = coeffs*(coeffs>0)
      negCoeffs = Array{Float64,1}(undef,length(indices))
      @. negCoeffs = coeffs*(coeffs<0)

      test = varLoBs[indices]
      test = varUpBs[indices]

      # compute the max and the min value of the constraint components
      maxArray = Array{Float64,1}(undef,length(indices))
      @. maxArray = posCoeffs*varUpBs[indices] + negCoeffs*varLoBs[indices]
      minArray = Array{Float64,1}(undef,length(indices))
      @. minArray = posCoeffs*varLoBs[indices] + negCoeffs*varUpBs[indices]


      # column index
      iteration = 0
      lastChanged = 1
      while true

            # select the variable to work on
            i = mod(iteration,length(indices))+1

            # compute new lower and upper bounds
            if coeffs[i] > 0
                  newLoB = (-(sum(maxArray[1:i-1]) + sum(maxArray[i+1:end])) + cnsLoB)/coeffs[i]
                  newUpB = (-(sum(minArray[1:i-1]) + sum(minArray[i+1:end])) + cnsUpB)/coeffs[i]
            else
                  newLoB = (-(sum(maxArray[1:i-1]) + sum(maxArray[i+1:end])) + cnsLoB)/coeffs[i]
                  newUpB = (-(sum(minArray[1:i-1]) + sum(minArray[i+1:end])) + cnsUpB)/coeffs[i]
            end

            if newLoB > newUpB
                  # Found infeasibility!
                  @info varLobs, varUpBs, newLoB, newUpB
                  error("Infeasible")
            end

            # perform changes
            if varLoBs[indices[i]] < newLoB
                  if indices[i] in dscIndices
                      newLoB = ceil(newLoB)
                  end

                  # change max and min arrays
                  if coeffs[i] > 0
                        minArray[i] = coeffs[i]*newLoB
                  else
                        maxArray[i] = coeffs[i]*newLoB
                  end
                  # update the bound
                  varLoBs[indices[i]] = newLoB
                  # remember that there was an update
                  lastChanged = i
                  # collect the updated variables
                  push!(updatedVars,indices[i])
            end

            if varUpBs[indices[i]] > newUpB
                  if indices[i] in dscIndices
                      newLoB = floor(newUpB)
                  end

                  # change max and min arrays
                  if coeffs[i] > 0
                        maxArray[i] = coeffs[i]*newUpB
                  else
                        minArray[i] = coeffs[i]*newUpB
                  end
                  # update the bound
                  varUpBs[indices[i]] = newUpB
                  # remember that there was an update
                  lastChanged = i
                  # collect the updated variables
                  push!(updatedVars,indices[i])
            end

            # If a variable bound changed, change cns bounds
            if lastChanged == i
                newCnsUpB = sum(maxArray[1:end])
                newCnsLoB = sum(minArray[1:end])
                if newCnsUpB < cnsUpB
                    cnsUpB = newCnsUpB
                end
                if newCnsLoB > cnsLoB
                    cnsLoB = newCnsLoB
                end
            end

            # if an update took place:
            # mark the current variable as updated and restart the iteration
            if lastChanged - 1 == mod(i,length(indices))
                  return updatedVars
            else
                  iteration = iteration + 1
            end
      end
end


# multi-constraint linear bounds propagation
function bounds_propagation!(rowsToCheck::Set{Int},
                            A::SparseMatrixCSC{Float64,Int},
                            cnsLoBs::Array{Float64,1},cnsUpBs::Array{Float64,1},
                            varLoBs::Array{Float64,1},varUpBs::Array{Float64,1},
                            dscIndices::Array{Int64,1})::Set{Int}

      # use a set for the output
      updatedVars = Set{Int}()

      while length(rowsToCheck) > 0

            # pick the first variable to update
            row = pop!(rowsToCheck)

            # perform bound propagation and collect the updated variables
            newUpdatedVars = bounds_propagation!(row, A, cnsLoBs, cnsUpBs, varLoBs, varUpBs, dscIndices)

            if length(newUpdatedVars) > 0
                  # collect the new rows to check
                  newRowsToCheck = unique(findnz(A[1:end,collect(newUpdatedVars)])[1])
                  deleteat!(newRowsToCheck,findfirst(x->x==row,newRowsToCheck)) # do not reinsert the current row
                  union!(rowsToCheck,newRowsToCheck)
                  # remember which variables were updated
                  union!(updatedVars,newUpdatedVars)
            end
      end

      return updatedVars
end

function bounds_propagation!(A::SparseMatrixCSC{Float64,Int},
                            cnsLoBs::Array{Float64,1},cnsUpBs::Array{Float64,1},
                            varLoBs::Array{Float64,1},varUpBs::Array{Float64,1},
                            dscIndices::Array{Int64,1})::Set{Int}
    return bounds_propagation!(Set(1:size(cnsLoBs)[1]),
                               A, cnsLoBs, cnsUpBs,
                               varLoBs, varUpBs, dscIndices)
end

function variable_bounds_propagation!(updatedVars::Array{Int64,1},
                            A::SparseMatrixCSC{Float64,Int},
                            cnsLoBs::Array{Float64,1},cnsUpBs::Array{Float64,1},
                            varLoBs::Array{Float64,1},varUpBs::Array{Float64,1},
                            dscIndices::Array{Int64,1})::Set{Int}
    return bounds_propagation!(Set(unique(findnz(A[1:end,collect(updatedVars)]))[1]),
                               A, cnsLoBs, cnsUpBs,
                               varLoBs, varUpBs, dscIndices)
end
