# @Author: Massimo De Mauri <massimo>
# @Date:   2019-06-19T16:35:26+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: python_interface.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-20T15:57:55+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# import the needed components
from os import path
from copy import copy, deepcopy
from warnings import warn
from numpy import array, matrix
from julia import Julia,OpenBB

class OpenBBmodel:

    def __init__(self):
        self.jl = OpenBB
        # current_directory = path.abspath(path.dirname(__file__))
        # self.jl.eval("using OpenBB")

    ######################## setup ########################
    # setup a OpenBB workspace
    def setup(self,problemDict=None,options={}):

        # get the subsolver name
        opts = copy(options)
        if not "subsolver" in opts:
            raise NameError("OpenBB, no subsolver name provided")
        else:
            subsolver = opts.pop("subsolver")

        # select a subsolver and load the list of available settings
        bb_settings_list = self.jl.eval_string("fieldnames(BBsettings)")
        if subsolver == "osqp":
            ss_settings_list = self.jl.eval_string("fieldnames(OSQPsettings)")
        elif subsolver == "gurobi":
            ss_settings_list = self.eval_string("fieldnames(GUROBIsettings)")
        elif subsolver == "qpalm":
            ss_settings_list = self.jl.eval_string("fieldnames(QPALMsettings)")
        else:
            raise NameError("OpenBB, unknown subsolver: "+subsolver)

        # create the settings dictionaries
        bb_settings = {}
        ss_settings = {}
        for k in opts:
            if k[:len(subsolver)+1] == subsolver+".":
                if k[len(subsolver)+1:] in ss_settings_list:
                    ss_settings[k[len(subsolver)+1:]] = opts[k]
                else:
                    warn("OpenBB, subsolver option not recognized: "+k[len(subsolver)+1:])
            else:
                if k in bb_settings_list:
                    bb_settings[k] = opts[k]
                else:
                    warn("OpenBB, option not recognized: "+k)

        # finally setup the current_workspace
        if problemDict is None:
            # set the opts for branch and bound
            self.jl.setup(subsolver,{},bb_settings,ss_settings)
        else:
            # check the input
            if not "objFun" in problemDict:
                raise NameError("keyword \"objFun\" not found")
            elif not "cnsSet" in problemDict:
                raise NameError("keyword \"cnsSet\" not found")
            elif not "varSet" in problemDict:
                raise NameError("keyword \"varSet\" not found")
            else:

                # reformat the objective
                localProblemDict = deepcopy(problemDict)
                if "Q" in localProblemDict["objFun"]: localProblemDict["objFun"]["Q"] = matrix(localProblemDict["objFun"]["Q"])
                localProblemDict["objFun"]["L"] = array(localProblemDict["objFun"]["L"]).flatten()

                # reformat the constraints set
                localProblemDict["cnsSet"]["A"] = matrix(localProblemDict["cnsSet"]["A"])
                localProblemDict["cnsSet"]["loBs"] = array(localProblemDict["cnsSet"]["loBs"]).flatten()
                localProblemDict["cnsSet"]["upBs"] = array(localProblemDict["cnsSet"]["upBs"]).flatten()

                # reformat the variables set
                localProblemDict["varSet"]["loBs"] = array(localProblemDict["varSet"]["loBs"]).flatten()
                localProblemDict["varSet"]["upBs"] = array(localProblemDict["varSet"]["upBs"]).flatten()
                if "vals" in localProblemDict["varSet"]: localProblemDict["varSet"]["vals"] = array(localProblemDict["varSet"]["vals"]).flatten()
                if "dscIndices" in localProblemDict["varSet"]: localProblemDict["varSet"]["dscIndices"] = [int(localProblemDict["varSet"]["dscIndices"][k])+1 for k in range(len(localProblemDict["varSet"]["dscIndices"]))]
                if "sos1Groups" in localProblemDict["varSet"]: localProblemDict["varSet"]["sos1Groups"] = [int(localProblemDict["varSet"]["sos1Groups"][k]) for k in range(len(localProblemDict["varSet"]["sos1Groups"]))]
                if "pseudoCosts" in localProblemDict["varSet"]: localProblemDict["varSet"]["pseudoCosts"] = array(localProblemDict["varSet"]["pseudoCosts"]).flatten()


            # create a branch and bound current_workspace
            self.jl.setup(subsolver,localProblemDict,bb_settings,ss_settings)
            return


    ######################## solve ########################
    # solve the problem
    def solve(self):
        # solve the problem
        self.jl.solve_b()
        return



    ######################## inspect ########################

    # ...
    def get_all(self):
        return self.jl.globalWorkspace

    # ...
    def get_settings(self):
        return self.jl.get_settings()

    # ...
    def get_subsolver_settings(self):
        return self.jl.get_subsolver_settings()

    # get statistics and results of the last optimization
    def get_status(self):
        return self.jl.get_status()

    # ...
    def get_best_solution(self,localOnly=False):
        return self.jl.get_best_solution(localOnly)

    # ...
    def get_all_solutions(self,localOnly=False):
        return self.jl.get_all_solutions(localOnly)

    # ...
    def get_best_node(self,localOnly=False):
        return self.jl.get_best_node(localOnly)

    # ...
    def get_numVariables(self):
        return self.jl.get_numVariables()

    # ...
    def get_numConstraints(self):
        return self.jl.get_numConstraints()

    # ...
    def get_numDiscreteVariables(self):
        return self.jl.get_numDiscreteVariables()


    # ...
    def get_constraints(self):
        return self.jl.get_constraints()

    # ...
    def get_objective(self):
        return self.jl.get_objective()

    # ...
    def get_constraints_sparsity(self):
        out = self.jl.get_constraints_sparsity()
        for k in range(len(out[0])):
            out[0][k] = out[0][k] - 1
            out[1][k] = out[1][k] - 1
        return out

    # ...
    def get_constraint_sparsity(self,index):
        out = self.jl.get_constraint_sparsity(index+1)
        for k in range(len(out)):
            out[k] = out[k] - 1
        return out

    # ...
    def get_objective_sparsity(self):
        out = self.jl.get_objective_sparsity()
        for k in range(len(out[0])):
            out[0][k] = out[0][k] - 1
            out[1][k] = out[1][k] - 1
        return out

    # ...
    def get_variableBounds(self):
        return self.jl.get_variableBounds()

    # ...
    def get_constraintBounds(self):
        return self.jl.get_constraintBounds()

    # ...
    def get_numActiveNodes(self):
        return self.jl.get_numActiveNodes()

    # ...
    def get_numUnactiveNodes(self):
        return self.jl.get_numUnactiveNodes()

    # ...
    def get_numSolutions(self):
        return self.jl.get_numSolutions()

    ######################## update workspace ########################
    # ...
    def reset_explored_nodes(self,localOnly=False):
        self.jl.reset_explored_nodes_b(localOnly)
        return

    # ...
    def update(self,localOnly=False):
        self.jl.update_b(localOnly)
        return

    # ...
    def reset(self,localOnly=False):
        self.jl.reset_b(localOnly)
        return

    # ...
    def clear(self,localOnly=False):
        self.jl.clear_b(localOnly)
        return


    ######################## update problem ########################
    # ...
    def append_constraints(self,constraintsDict,suppressWarnings=False,suppressUpdate=False,localOnly=False):
        if not "A" in constraintsDict or not "loBs" in constraintsDict or not "upBs" in constraintsDict:
            raise NameError("constraintsDict has to be a dictionary with the following keywords: A,loBs,upBs")
        self.jl.append_constraints_b(constraintsDict,suppressWarnings,suppressUpdate,localOnly)
        return

    # ...
    def insert_constraints(self,constraintsDict,index,suppressWarnings=False,suppressUpdate=False,localOnly=False):
        if not "A" in constraintsDict or not "loBs" in constraintsDict or not "upBs" in constraintsDict:
            raise NameError("constraintsDict has to be a dictionary with the following keywords: A,loBs,upBs")
        self.jl.insert_constraints_b(constraintsDict,index+1,suppressWarnings,suppressUpdate,localOnly)
        return


    # ...
    def remove_constraints(self,indices,suppressWarnings=False,suppressUpdate=False,localOnly=False):
        for k in range(len(indices)):
            indices[k] = indices[k] + 1
        self.jl.remove_constraints_b(indices,suppressWarnings,suppressUpdate,localOnly)
        return

    # ...
    def permute_constraints(self,permutation,suppressWarnings=False,suppressUpdate=False,localOnly=False):
        for k in range(len(permutation)):
            permutation[k] = permutation[k] + 1
        self.jl.permute_constraints_b(permutation,suppressWarnings,suppressUpdate,localOnly)
        return

    # ...
    def update_bounds(self,boundsDict,suppressWarnings=False,suppressUpdate=False,localOnly=False):
        self.jl.update_bounds_b(boundsDict,suppressWarnings,suppressUpdate,localOnly)
        return


    # adds new variables, new constraints and a new part of the objective to the problem
    def append_problem(self,problemDict,suppressWarnings=False,suppressUpdate=False,localOnly=False):
        if not problemDict is None:
            # check the input
            if not "objFun" in problemDict:
                raise NameError("keyword \"objFun\" not found")
            elif not "cnsSet" in problemDict:
                raise NameError("keyword \"cnsSet\" not found")
            elif not "varSet" in problemDict:
                raise NameError("keyword \"varSet\" not found")
            else:

                # reformat the objective
                localProblemDict = deepcopy(problemDict)
                if "Q" in localProblemDict["objFun"]: localProblemDict["objFun"]["Q"] = matrix(localProblemDict["objFun"]["Q"])
                localProblemDict["objFun"]["L"] = array(localProblemDict["objFun"]["L"]).flatten()

                # reformat the constraints set
                localProblemDict["cnsSet"]["A"] = matrix(localProblemDict["cnsSet"]["A"])
                localProblemDict["cnsSet"]["loBs"] = array(localProblemDict["cnsSet"]["loBs"]).flatten()
                localProblemDict["cnsSet"]["upBs"] = array(localProblemDict["cnsSet"]["upBs"]).flatten()

                # reformat the variables set
                localProblemDict["varSet"]["loBs"] = array(localProblemDict["varSet"]["loBs"]).flatten()
                localProblemDict["varSet"]["upBs"] = array(localProblemDict["varSet"]["upBs"]).flatten()
                if "vals" in localProblemDict["varSet"]: localProblemDict["varSet"]["vals"] = array(localProblemDict["varSet"]["vals"]).flatten()
                if "dscIndices" in localProblemDict["varSet"]: localProblemDict["varSet"]["dscIndices"] = [int(localProblemDict["varSet"]["dscIndices"][k])+1 for k in range(len(localProblemDict["varSet"]["dscIndices"]))]
                if "sos1Groups" in localProblemDict["varSet"]: localProblemDict["varSet"]["sos1Groups"] = [int(localProblemDict["varSet"]["sos1Groups"][k]) for k in range(len(localProblemDict["varSet"]["sos1Groups"]))]
                if "pseudoCosts" in localProblemDict["varSet"]: localProblemDict["varSet"]["pseudoCosts"] = array(localProblemDict["varSet"]["pseudoCosts"]).flatten()

                self.jl.append_problem_b(localProblemDict,suppressWarnings,suppressUpdate,localOnly)
        return


    # marks as discrete formally not discrete variables
    def integralize_variables(self,newDscIndices,newSos1Groups=array([],int),suppressWarnings=False,suppressUpdate=False,localOnly=False):
        for k in range(len(newDscIndices)):
            newDscIndices[k] = newDscIndices[k] + 1
        self.jl.integralize_variables_b(newDscIndices,newSos1Groups,suppressWarnings,suppressUpdate,localOnly)
        return

    # ...
    def update_objectiveCutoff(self,newCutoff,suppressWarnings=False,suppressUpdate=False,localOnly=False):
        self.jl.update_objectiveCutoff_b(newCutoff,suppressWarnings,suppressUpdate,localOnly)
        return
