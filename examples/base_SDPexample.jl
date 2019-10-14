# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:22:14+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename:
# @Last modified by:   massimo
# @Last modified time: 2019-07-15T12:22:31+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

using OpenBB
using SparseArrays
using LinearAlgebra
using Juno

subsolver = "JUMP_SDP"
# subsolver = "gurobi"
if subsolver == "JUMP_SDP"
    subsolverSettings = OpenBB.JUMP_SDPsettings()
elseif subsolver == "gurobi"
     subsolverSettings = OpenBB.GUROBIsettings()
end
B = rand(4,4);
A = Symmetric(B);
g = 5*rand(4,4);
L = Symmetric(g);

# Basic usage of OpenBB for mixed-integer quadratic problems
problem = OpenBB.Problem(objFun=OpenBB.LinearObjective(L),
                         cnsSet=OpenBB.SDPConstraintSet(A),
                         varSet=OpenBB.VariableSet(loBs=[-Infs(4,4)],upBs=[Infs(4,4)],vals = zeros(4,4)))


workspace = OpenBB.setup(problem,OpenBB.BBsettings(verbose=true,interactiveMode=true,numProcesses=1),subsolverSettings)
#result = OpenBB.solve!(workspace)
