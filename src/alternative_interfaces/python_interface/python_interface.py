# @Author: Massimo De Mauri <massimo>
# @Date:   2019-06-19T16:35:26+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: python_interface.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-19T19:44:07+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# import the needed components
from os import path
from copy import copy
from warnings import warn
from numpy import array, matrix
from julia import Julia

class OpenBBmodel:

    def __init__(self):
        self.jl = Julia()
        current_directory = path.abspath(path.dirname(__file__))
        self.jl.eval("include(\""+current_directory+"/../ctypes_interface.jl\")")

    ######################## setup ########################
    # setup a OpenBB current_workspace
    def setup(self,problemDict=None,options={}):

        # get the subsolver name
        opts = copy(options)
        if not "subsolver" in opts:
            raise NameError("OpenBB, no subsolver name provided")
        else:
            subsolver = opts.pop("subsolver")

        # select a subsolver and load the list of available settings
        bb_settings_list = self.jl.fieldnames(self.jl.eval("OpenBB.BBsettings"))
        if subsolver == "osqp":
            ss_settings_list = self.jl.fieldnames(self.jl.eval("OpenBB.OSQPsettings"))
        elif subsolver == "gurobi":
            ss_settings_list = self.jl.fieldnames(self.jl.eval("OpenBB.GUROBIsettings"))
        elif subsolver == "qpalm":
            ss_settings_list = self.jl.fieldnames(self.jl.eval("OpenBB.QPALMsettings"))
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
                if "Q" in problemDict["objFun"]: problemDict["objFun"]["Q"] = matrix(problemDict["objFun"]["Q"])
                problemDict["objFun"]["L"] = array(problemDict["objFun"]["L"]).flatten()

                # reformat the constraints set
                problemDict["cnsSet"]["A"] = matrix(problemDict["cnsSet"]["A"])
                problemDict["cnsSet"]["loBs"] = array(problemDict["cnsSet"]["loBs"]).flatten()
                problemDict["cnsSet"]["upBs"] = array(problemDict["cnsSet"]["upBs"]).flatten()

                # reformat the variables set
                problemDict["varSet"]["loBs"] = array(problemDict["varSet"]["loBs"]).flatten()
                problemDict["varSet"]["upBs"] = array(problemDict["varSet"]["upBs"]).flatten()
                if "vals" in problemDict["varSet"]: problemDict["varSet"]["vals"] = array(problemDict["varSet"]["vals"]).flatten()
                if "dscIndices" in problemDict["varSet"]: problemDict["varSet"]["dscIndices"] = [int(problemDict["varSet"]["dscIndices"][k])+1 for k in range(len(problemDict["varSet"]["dscIndices"]))]
                if "sos1Groups" in problemDict["varSet"]: problemDict["varSet"]["sos1Groups"] = [int(problemDict["varSet"]["sos1Groups"][k]) for k in range(len(problemDict["varSet"]["sos1Groups"]))]
                if "pseudoCosts" in problemDict["varSet"]: problemDict["varSet"]["pseudoCosts"] = array(problemDict["varSet"]["pseudoCosts"]).flatten()


            # create a branch and bound current_workspace
            self.jl.setup(subsolver,problemDict,bb_settings,ss_settings)
            return


    ######################## solve ########################
    # solve the problem
    def solve(self):
        # solve the problem
        self.jl.solve_b()
        return



    ######################## inspect ########################
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
    def get_constraints_sparsity(self):
        out = self.jl.get_constraints_sparsity()
        for k in range(len(out[1])):
            out[1][k] = out[1][k] - 1
            out[2][k] = out[2][k] - 1
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
        for k in range(len(out[1])):
            out[1][k] = out[1][k] - 1
            out[2][k] = out[2][k] - 1
        return out

    # ...
    def get_variableBounds(self):
        return self.jl.get_variableBounds()

    # ...
    def get_constraintBounds(self):
        return self.jl.get_constraintBounds()



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
        self.jl.remove_constraints_b(permutation,suppressWarnings,suppressUpdate,localOnly)
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
                if "Q" in problemDict["objFun"]: problemDict["objFun"]["Q"] = matrix(problemDict["objFun"]["Q"])
                problemDict["objFun"]["L"] = array(problemDict["objFun"]["L"]).flatten()

                # reformat the constraints set
                problemDict["cnsSet"]["A"] = matrix(problemDict["cnsSet"]["A"])
                problemDict["cnsSet"]["loBs"] = array(problemDict["cnsSet"]["loBs"]).flatten()
                problemDict["cnsSet"]["upBs"] = array(problemDict["cnsSet"]["upBs"]).flatten()

                # reformat the variables set
                problemDict["varSet"]["loBs"] = array(problemDict["varSet"]["loBs"]).flatten()
                problemDict["varSet"]["upBs"] = array(problemDict["varSet"]["upBs"]).flatten()
                if "vals" in problemDict["varSet"]: problemDict["varSet"]["vals"] = array(problemDict["varSet"]["vals"]).flatten()
                if "dscIndices" in problemDict["varSet"]: problemDict["varSet"]["dscIndices"] = [int(problemDict["varSet"]["dscIndices"][k])+1 for k in range(len(problemDict["varSet"]["dscIndices"]))]
                if "sos1Groups" in problemDict["varSet"]: problemDict["varSet"]["sos1Groups"] = [int(problemDict["varSet"]["sos1Groups"][k]) for k in range(len(problemDict["varSet"]["sos1Groups"]))]
                if "pseudoCosts" in problemDict["varSet"]: problemDict["varSet"]["pseudoCosts"] = array(problemDict["varSet"]["pseudoCosts"]).flatten()

                self.jl.append_problem_b(problemDict,suppressWarnings,suppressUpdate,localOnly)
        return


    # marks as discrete formally not discrete variables
    def integralize_variables(self,newDscIndices,suppressWarnings=False,suppressUpdate=False,localOnly=False):
        for k in range(len(newDscIndices)):
            newDscIndices[k] = newDscIndices[k] + 1
        self.jl.integralize_variables_b(newDscIndices,suppressWarnings,suppressUpdate,localOnly)
        return

    # ...
    def update_cutoff(self,newCutoff,suppressWarnings=False,suppressUpdate=False,localOnly=False):
        self.jl.update_cutoff_b(newCutoff,suppressWarnings,suppressUpdate,localOnly)
        return
