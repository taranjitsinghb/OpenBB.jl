# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-26T15:03:06+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: ConstraintSets.jl
# @Last modified by:   massimo
# @Last modified time: 2019-05-23T18:56:26+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# abstract and null types
abstract type AbstractConstraintSet end


struct NullConstraintSet <: AbstractConstraintSet end

function get_numConstraints(constraintSet::NullConstraintSet)::Int
    return 0
end

import SparseArrays.sparse
function sparse(constraintSet::NullConstraintSet)::LinearConstraintSet
    return constraintSet
end

function get_sparsity(constraintSet::NullConstraintSet)::Tuple{Array{Int,1},Array{Int,1}}
    return (Int[],Int[])
end

include("./LinearConstraintSet.jl")
