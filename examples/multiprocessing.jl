# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-09T18:02:42+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: multiprocessing.jl
# @Last modified by:   massimo
# @Last modified time: 2019-05-23T19:17:13+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}



clearconsole()
using Distributed
using SharedArrays
using LinearAlgebra
using SparseArrays
using OpenBB






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
                         varSet=OpenBB.VariableSet(loBs=[-5.;-Infs(3)],upBs=[ 5.;Infs(3)],val=zeros(4),dscIndices=[1]))


workspace = OpenBB.setup(problem,OpenBB.BBsettings(verbose=true,dynamicMode=true),subsolverSettings)
# OpenBB.solve!(workspace)

#
# for p in workers()
# 	@spawnat p problem,OpenBB.BBsettings(verbose=true,dynamicMode=true),subsolverSettings
# end
#
# @everywhere workspace = OpenBB.setup(problem,OpenBB.BBsettings(verbose=true,dynamicMode=true),subsolverSettings)
#
# @everywhere OpenBB.run!(workspace)
