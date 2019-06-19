# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-16T15:53:30+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: OSQP_inspect.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-19T21:13:27+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}



# this function returns the sparsity pattern of the constraint set
function get_constraints_sparsity(workspace::OSQPworkspace)::Tuple{Array{Int,1},Array{Int,1}}
    return findnz(sparse(workspace.A))[1:2]
end

# this function returns the sparsity pattern a constraint in the constraint set
function get_constraint_sparsity(workspace::OSQPworkspace,index::Int)::Array{Int,1}
    return findnz(sparse(workspace.A[index,:]))[1]
end

# this function returns the sparsity pattern of the objevtive function
function get_objective_sparsity(workspace::OSQPworkspace)::Tuple{Array{Int,1},Array{Int,1}}
    return findnz(sparse(workspace.Q))[1:2]
end

#
function get_numVariables(workspace::OSQPworkspace)::Int
    return size(workspace.A,2)
end

#
function get_numConstraints(workspace::OSQPworkspace)::Int
    return size(workspace.A,1)
end


#
function get_variableBounds(workspace::OSQPworkspace)::Tuple{Array{Float64,1},Array{Float64,1}}
    return (workspace.varLoBs,workspace.varUpBs)
end

#
function get_constraintBounds(workspace::OSQPworkspace)::Tuple{Array{Float64,1},Array{Float64,1}}
    return (workspace.cnsLoBs,workspace.cnsUpBs)
end

# ...
function get_settings(workspace::OSQPworkspace)::OSQPsettings
    return workspace.settings
end
