# @Author: Massimo De Mauri <massimo>
# @Date:   2019-08-27T12:38:57+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: test_problem_definition_update.jl
# @Last modified by:   massimo
# @Last modified time: 2019-08-27T12:39:24+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

using OpenBB
using LinearAlgebra
using SparseArrays


# test variable set
for v in var
    tmp = deepcopy(v)
    OpenBB.insert_variables!(tmp,tmp,4)
    @assert OpenBB.get_size(tmp) == 8
    OpenBB.remove_variables!(tmp,[2,4])
    @assert OpenBB.get_size(tmp) == 6
    @assert all(@. tmp.dscIndices == [1,2,4])
    @assert length(tmp.dscIndices) == length(tmp.sos1Groups) == size(tmp.pseudoCosts[1],1) == size(tmp.pseudoCosts[2],1)
    @assert all(@. tmp.sos1Groups == [0,1,2])
end

# test constraint sets
for c in cns
    tmp = deepcopy(c)
    OpenBB.update_bounds!(tmp,[2.],[2.])
    OpenBB.insert_constraints!(tmp,c,2)
    OpenBB.append_constraints!(tmp,tmp)
    OpenBB.remove_constraints!(tmp,[3])
    OpenBB.append_variables!(tmp,5)
    OpenBB.remove_variables!(tmp,[1])
    OpenBB.permute_constraints!(tmp,[2,1,3])

    @assert OpenBB.get_numVariables(tmp) == 8
    @assert OpenBB.get_size(tmp) == 3
    @assert OpenBB.get_sparsity(tmp) == ([1, 2, 3], [1, 1, 1])
    @assert OpenBB.get_bounds(tmp) == ([1.0, 2.0, 1.0], [1.0, 2.0, 1.0])
end




# test problem setup
for o in obj
    for c in cns
        for v in var
            tmp = OpenBB.Problem(objFun=o,cnsSet=c,varSet=v)
        end
    end
end

println(" - problem definition update, ok")
