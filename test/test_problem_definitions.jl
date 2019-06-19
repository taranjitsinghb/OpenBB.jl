# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-05T15:04:00+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: test_problem_definitions.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-20T00:49:22+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

obj = Array{OpenBB.AbstractObjectiveFunction,1}(undef,3)
obj[1] = OpenBB.QuadraticObjective(Q=Matrix(1.0I,4,4,),L=ones(4))
obj[2] = OpenBB.LinearObjective(L=ones(4))
obj[3] = OpenBB.NullObjectiveFunction()


cns = Array{OpenBB.AbstractConstraintSet,1}(undef,2)
cns[1] = OpenBB.LinearConstraintSet(A=ones(1,4),loBs=[1.],upBs=[1.])
cns[2] = OpenBB.NullConstraintSet()

var = Array{OpenBB.VariableSet,1}(undef,2)
var[1] = OpenBB.VariableSet(loBs=zeros(4),upBs=Infs(4),vals=zeros(4),dscIndices=[1])
var[2] = OpenBB.EmptyVarSet()

for o in obj
    for c in cns
        for v in var
             OpenBB.Problem(objFun=o,cnsSet=c,varSet=v)
        end
    end
end

println(" - problem definitions, ok")
