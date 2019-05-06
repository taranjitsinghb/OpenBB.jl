# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-05T15:04:00+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: test_problem_definitions.jl
# @Last modified by:   massimo
# @Last modified time: 2019-03-05T18:44:17+01:00
# @License: apache 2.0
# @Copyright: {{copyright}}

obj = Array{OpenBB.AbstractObjectiveFun,1}(undef,3)
obj[1] = OpenBB.QuadraticObj(Q=Matrix(1.0I,4,4,),L=ones(4))
obj[2] = OpenBB.LinearObj(L=ones(4))
obj[3] = OpenBB.NullObjectiveFun()


cns = Array{OpenBB.AbstractConstraintSet,1}(undef,2)
cns[1] = OpenBB.LinearCns(A=ones(1,4),loBs=[1.],upBs=[1.],sosIndices=[1])
cns[2] = OpenBB.NullConstraintSet()

var = Array{OpenBB.VariableSet,1}(undef,2)
var[1] = OpenBB.VariableSet(loBs=zeros(4),upBs=Infs(4),val=zeros(4),dscIndices=[1])
var[2] = OpenBB.EmptyVarSet()

for o in obj
    for c in cns
        for v in var
             OpenBB.Problem(objFun=o,cnsSet=c,varSet=v)
        end
    end
end

println(" - problem definitions, ok")
