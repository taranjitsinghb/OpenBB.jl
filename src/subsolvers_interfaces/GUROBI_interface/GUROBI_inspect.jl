# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-16T15:47:57+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: GUROBI_inspect.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-20T15:15:29+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}



# this function returns the sparsity pattern of the constraint set
function get_constraints_sparsity(workspace::GUROBIworkspace)::Tuple{Array{Int,1},Array{Int,1}}
    return findnz(sparse(workspace.A))[1:2]
end

# this function returns the sparsity pattern a constraint in the constraint set
function get_constraint_sparsity(workspace::GUROBIworkspace,index::Int)::Array{Int,1}
    return findnz(sparse(workspace.A[index,:]))[1]
end


# this function returns the sparsity pattern of the objective function
function get_objective_sparsity(workspace::GUROBIworkspace)::Tuple{Array{Int,1},Array{Int,1}}
    return findnz(sparse(workspace.Q))[1:2]
end

#
function get_numVariables(workspace::GUROBIworkspace)::Int
    return size(workspace.A,2)
end

#
function get_numConstraints(workspace::GUROBIworkspace)::Int
    return size(workspace.A,1)
end


#
function get_variableBounds(workspace::GUROBIworkspace)::Tuple{Array{Float64,1},Array{Float64,1}}
    return (workspace.varLoBs,workspace.varUpBs)
end

#
function get_constraintBounds(workspace::GUROBIworkspace)::Tuple{Array{Float64,1},Array{Float64,1}}
    return (workspace.cnsLoBs,workspace.cnsUpBs)
end

# ...
function get_settings(workspace::GUROBIworkspace)::GUROBIsettings
    return workspace.settings
end

# ...
function get_constraints(workspace::GUROBIworkspace)::LinearConstraintSet
    return LinearConstraintSet(A=workspace.A,loBs=workspace.cnsLoBs,upBs=workspace.cnsUpBs)
end

# ...
function get_objective(workspace::GUROBIworkspace)::QuadraticObjective
    return QuadraticObjective(Q=workspace.Q,L=workspace.L)
end
