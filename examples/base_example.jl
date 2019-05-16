# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:22:14+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename:
# @Last modified by:   massimo
# @Last modified time: 2019-05-15T17:33:14+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}

using OpenBB
using SparseArrays
using LinearAlgebra

clearconsole()

subsolver = "osqp"
# subsolver = "gurobi"
if subsolver == "osqp"
    subsolverSettings = OpenBB.OSQPsettings()
elseif subsolver == "gurobi"
     subsolverSettings = OpenBB.GUROBIsettings()
end


# Basic usage of OpenBB for mixed-integer quadratic problems
problem = OpenBB.Problem(objFun=OpenBB.QuadraticObj(Q=Matrix(2.0I,4,4,),L=[-.5,0.,0.,0.]),
                         cnsSet=OpenBB.LinearCns(A=ones(1,4),loBs=[1.],upBs=[1.]),
                         varSet=OpenBB.VariableSet(loBs=[-5.;-Infs(3)],upBs=[ 5.;Infs(3)],val=zeros(4),dscIndices=[1]))


workspace = OpenBB.setup(problem,OpenBB.BBsettings(verbose=true,dynamicMode=true,numProcesses=2),subsolverSettings)
result = OpenBB.solve!(workspace)
