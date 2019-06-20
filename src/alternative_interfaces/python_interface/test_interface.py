# @Author: Massimo De Mauri <massimo>
# @Date:   2019-06-19T18:15:50+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: test.py
# @Last modified by:   massimo
# @Last modified time: 2019-06-20T16:00:33+02:00
# @License: LGPL-3.0

import numpy as np
from python_interface import *


# load openBB
workspace = OpenBBmodel()


# create the problem
problem = {"objFun":{},"cnsSet":{},"varSet":{}}

problem["objFun"]["Q"] = 2*np.eye(5)
problem["objFun"]["Q"][2,2] = 0
problem["objFun"]["L"] = np.zeros(5)

problem["cnsSet"]["A"] = np.matrix([[1.,1.,0.,-1.,-1.],[0.,1.,1.,1.,0.]])
problem["cnsSet"]["loBs"] = np.array([0.,1.])
problem["cnsSet"]["upBs"] = np.array([0.,1.])

problem["varSet"]["loBs"] = [.5,0,0,0,-10]
problem["varSet"]["upBs"] = [.5,1,1,1,10]
problem["varSet"]["vals"] = [.5,0,1,0,.5]
problem["varSet"]["dscIndices"] = [1,2,3]
problem["varSet"]["sos1Groups"] = [1,1,1]

settings = {"subsolver":"osqp",
            "verbose":False,
            "numProcesses":1,
            "dynamicMode":True
            }

workspace.setup(problem,settings)
workspace.solve()

assert workspace.get_best_solution()["objective"] >= 0.5 - workspace.get_settings()["primalTolerance"]
assert workspace.get_best_solution()["objective"] <= 0.5 + workspace.get_settings()["primalTolerance"]

workspace.get_all_solutions()
workspace.get_best_node()

assert workspace.get_numVariables() == 5
assert workspace.get_numConstraints() == 2
assert workspace.get_numDiscreteVariables() == 3

workspace.get_constraints_sparsity()
assert all(workspace.get_constraint_sparsity(0) == [0,1,3,4])
assert all(workspace.get_constraint_sparsity(1) == [1,2,3])
assert all(workspace.get_objective_sparsity()[0] == [0, 1, 3, 4])

assert all(workspace.get_variableBounds()[0] == [0.5, 0.0, 0.0, 0.0, -10.0])
assert all(workspace.get_constraintBounds()[0] == [0.0, 1.0])
assert workspace.get_status()["description"] == "optimalSolutionFound"

workspace.reset_explored_nodes()
assert workspace.get_numUnactiveNodes() == workspace.get_numSolutions() == 0
assert workspace.get_status()["objUpB"] == np.inf

workspace.clear()
assert workspace.get_numActiveNodes() == 0

workspace.reset()
assert workspace.get_numActiveNodes() == 1

workspace.solve()
workspace.append_constraints(problem["cnsSet"],False,True,False)
workspace.remove_constraints([2,3],True,True,False)
workspace.permute_constraints([1,0],True,True,False)
assert all(workspace.get_constraintBounds()[0] == [1.0, 0.0])


workspace.update_bounds({"cnsLoBs":[0.0,0.0],"varLoBs":[0.5, 0.0, 0.0, 0.0, -100.0]},False,True,False)
assert all(workspace.get_constraintBounds()[0] == [0.0, 0.0])
assert all(workspace.get_variableBounds()[0] == [0.5, 0.0, 0.0, 0.0, -100.0])
workspace.update_bounds({"cnsLoBs":[1.0,0.0]},False,True,False)


workspace.append_problem(problem,False,False,False)
workspace.solve()
assert workspace.get_best_solution()["objective"] >= 1.0 - workspace.get_settings()["primalTolerance"]
assert workspace.get_best_solution()["objective"] <= 1.0 + workspace.get_settings()["primalTolerance"]


workspace.integralize_variables([4])
workspace.solve()
assert workspace.get_status()["description"] == "infeasible"


workspace.get_constraints()
workspace.get_objective()

print(" - Python interface, ok")
