# @Author: Massimo De Mauri <massimo>
# @Date:   2019-06-19T18:15:50+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: test.py
# @Last modified by:   massimo
# @Last modified time: 2019-06-20T00:42:15+02:00
# @License: LGPL-3.0

import numpy as np
from python_interface import *


# load openBB
workspace = OpenBBmodel()


# create the problem
problem = {"objFun":{},"cnsSet":{},"varSet":{}}

problem["objFun"]["Q"] = .5*np.eye(5)
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
            "verbose":True,
            "numProcesses":2,
              }

workspace.setup(problem,settings)
workspace.solve()
