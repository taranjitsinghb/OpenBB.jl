# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-05T15:01:45+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: test_QP.jl
# @Last modified by:   massimo
# @Last modified time: 2019-07-15T12:43:53+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

function test_QP_subsolver(subsolver)

    if subsolver == "osqp"
        subsolverSettings = OpenBB.OSQPsettings()
    elseif subsolver == "gurobi"
        subsolverSettings = OpenBB.GUROBIsettings()
    elseif subsolver == "qpalm"
        subsolverSettings = OpenBB.QPALMsettings()
    end
    print(" - ")
    print("1...")
    # create first problem
    problem = OpenBB.Problem(objFun=OpenBB.QuadraticObjective(Q=Matrix(1.0I,4,4,),L=[-.5,0.,0.,0.]),
                             cnsSet=OpenBB.LinearConstraintSet(A=ones(0,4),loBs=Float64[],upBs=Float64[]),
                             varSet=OpenBB.VariableSet(loBs=[-5.;-Infs(3)],upBs=[ 5.;Infs(3)],vals=zeros(4),dscIndices=[1]))
    workspace = OpenBB.setup(problem,OpenBB.BBsettings(interactiveMode=true,verbose=false,statusInfoPeriod=0.01,numProcesses=1),subsolverSettings)
    result0 = OpenBB.solve!(workspace)


    print("2...")
    # add some linear contraints
    OpenBB.append_constraints!(workspace,OpenBB.LinearConstraintSet(ones(1,4),[1.],[1.]))
    result1 = OpenBB.solve!(workspace)

    print("3...")
    # Basic usage of OpenBB for mixed-integer quadratic problems
    problem2 = OpenBB.Problem(objFun=OpenBB.QuadraticObjective(Q=sparse([1,2,3,4],[1,2,3,4],[1.,2.,3.,4.]),L=[2.,2.,2.,2.]),
                             cnsSet=OpenBB.LinearConstraintSet(A=ones(1,4),loBs=[1.],upBs=[1.]),
                             varSet=OpenBB.VariableSet(loBs=[-5.;-Infs(3)],upBs=[ 5.;Infs(3)],vals=zeros(4),dscIndices=[1]))

    OpenBB.append_problem!(workspace,problem2,suppressUpdate=true)
    OpenBB.permute_constraints!(workspace,reverse(collect(1:size(workspace.subsolverWS.A,1))),suppressUpdate=true)
    OpenBB.update!(workspace)
    result2 = OpenBB.solve!(workspace)

    print("4...")
    OpenBB.update_bounds!(workspace;varLoBs=[1.,0.,0.,0.,1.,0.,0.,0.])
    result3 = OpenBB.solve!(workspace)


    OpenBB.get_constraints_sparsity(workspace.subsolverWS)
    OpenBB.get_objective_sparsity(workspace.subsolverWS)

    println(subsolver,": setup + solve + update, ok")
end

if OpenBB.withOSQP()
    test_QP_subsolver("osqp")
end
if OpenBB.withGUROBI()
    test_QP_subsolver("gurobi")
end
if OpenBB.withQPALM()
    test_QP_subsolver("qpalm")
end
