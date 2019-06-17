# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-11T17:03:23+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: linear_bounds_propagation.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-14T11:45:59+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}


# single constraint linear bounds propagation
function bounds_propagation!(cns::SparseVector{Float64,Int},
                            cnsLoB::Float64,cnsUpB::Float64,
                            varLoBs::Array{Float64,1},varUpBs::Array{Float64,1})::Set{Int}

      # use a set for the output
      updatedVars = Set{Int}()

      # get the sparsity of cns
      indices,coeffs = findnz(cns)

      # distiguish between positive and negative coefficients
      posCoeffs = Array{Float64,1}(undef,length(indices))
      @. posCoeffs = coeffs*(coeffs>0)
      negCoeffs = Array{Float64,1}(undef,length(indices))
      @. negCoeffs = coeffs*(coeffs<0)

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
                  println(varLoBs,varUpBs)
                  stop
            end


            # perform changes
            if varLoBs[indices[i]] < newLoB
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
                            cnsLoBs::Array{Float64,1},cnsUpBs::Array{Float64},
                            varLoBs::Array{Float64,1},varUpBs::Array{Float64,1})

      # use a set for the output
      updatedVars = Set{Int}()

      while length(rowsToCheck) > 0

            # pick the first variable to update
            row = pop!(rowsToCheck)

            # perform bound propagation and collect the updated variables
            newUpdatedVars = bounds_propagation!(A[row,1:end],cnsLoBs[row],cnsUpBs[row],varLoBs,varUpBs)

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
