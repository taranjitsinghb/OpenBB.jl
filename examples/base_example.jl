# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:22:14+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename:
# @Last modified by:   massimo
# @Last modified time: 2019-06-19T16:11:34+02:00
# @License: LGPL-3.0
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
problem = OpenBB.Problem(objFun=OpenBB.QuadraticObjective(Q=Matrix(2.0I,4,4,),L=[-.5,0.,0.,0.]),
                         cnsSet=OpenBB.LinearConstraintSet(A=ones(1,4),loBs=[1.],upBs=[1.]),
                         varSet=OpenBB.VariableSet(loBs=[-5.;-Infs(3)],upBs=[ 5.;Infs(3)],vals=zeros(4),dscIndices=[1]))


workspace = OpenBB.setup(problem,OpenBB.BBsettings(verbose=true,dynamicMode=true,numProcesses=2),subsolverSettings)
result = OpenBB.solve!(workspace)
