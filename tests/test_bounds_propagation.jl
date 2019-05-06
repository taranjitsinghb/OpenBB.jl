# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-11T20:08:25+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: test_bounds_propagation.jl
# @Last modified by:   massimo
# @Last modified time: 2019-04-04T00:06:31+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}


using Revise
using OpenBB
using LinearAlgebra
using SparseArrays


# lets test a simple sos1 constraint
cns = SparseVector(4,[1,2,3,4],ones(4))
cnsUpB = 1.
cnsLoB = 1.


# set one variable to 1
varLoBs = vcat(1,zeros(3))
varUpBs = vcat(1, ones(3))
desiredVarLoBs = vcat(1,zeros(3))
desiredVarUpBs = vcat(1,zeros(3))
fixedVars = OpenBB.bounds_propagation!(cns,cnsLoB,cnsUpB,varLoBs,varUpBs)

for k in 2:4
    @assert k in fixedVars
end
@assert all(@. varLoBs == desiredVarLoBs)
@assert all(@. varUpBs == desiredVarUpBs)

# set one variable to 0
varLoBs = vcat(0,zeros(3))
varUpBs = vcat(0, ones(3))
desiredVarLoBs = zeros(4)
desiredVarUpBs = vcat(0,ones(3))

fixedVars = OpenBB.bounds_propagation!(cns,cnsLoB,cnsUpB,varLoBs,varUpBs)

@assert length(fixedVars) == 0
@assert all(@. varLoBs == desiredVarLoBs)
@assert all(@. varUpBs == desiredVarUpBs)



# lets test a couple of constraints
cns = sparse([1,1,1,2,2],[1,2,3,3,4],[1.,1.,1.,1.,1.])
cnsUpBs = ones(2)
cnsLoBs = ones(2)

# set the first variable to 1
varLoBs = vcat(1,zeros(3))
varUpBs = vcat(1, ones(3))
desiredVarLoBs = vcat(1,zeros(2),1)
desiredVarUpBs = vcat(1,zeros(2),1)

fixedVars = OpenBB.bounds_propagation!(Set([1]),cns,cnsLoBs,cnsUpBs,varLoBs,varUpBs)
for k in 2:4
    @assert k in fixedVars
end
@assert all(@. varLoBs == desiredVarLoBs)
@assert all(@. varUpBs == desiredVarUpBs)
