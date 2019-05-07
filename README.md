# OpenBB
An open and modular Branch and Bound framework written in Julia.
* version 0.1.0

## The Idea
Those are the driving ideas behind this project:
* It should be very easy to define and tackle different problem classes. Possibly, via the addition of new subsolvers.
* It must be possible to interrupt the Branch and Bound execution at any moment and restart it later.
* The state of the Branch & Bound should be transparent to the user and easy to manipulate.
* It must be possible to modify the problem (as much as it makes sense) without restarting the Branch and Bound execution.
* The code should be optimized as much as possible without harming the generality

## What is already here
* A basic sequential Branch and Bound implementation for Mixed-Integer Quadratic Problems with linear constraints.
* OSQP interface (https://osqp.org/).
* Gurobi (QP only) interface (http://www.gurobi.com/) based on Gurobi.jl (https://github.com/JuliaOpt/Gurobi.jl)

The code allows you to easily define your custom priority rules for the selection of the next subproblem to solve and the next variable to branch on. Moreover, it is very easy to define custom stopping criteria for Branch & Bound. Moreover, the whole code is subsolver-invariant and a new subsolver can be added by simply overloading the interface functions that the user can find in the "subsolvers_interfaces" folder. Finally, OpenBB provides a number of functions that allow safe constraints addition, problem expansion, bounds restrictions, etc..



### Disclaimer 1
The project is a early phase of development. Some features are still missing. Here is a temptative to do list:
* Add MPC functionalities (my main reason for this project)
* Parallelize the code
* Add qpOASES interface (QP solver with MPC in mind: https://projects.coin-or.org/qpOASES)
* Add SDP support
* Add NLP support (maybe... as usually Outer Approximation outperforms NLP B&B)
* Add bound propagation and cuts generation techniques
* Add presolvers
* Document all

Feel free to propose missing functionalities and to collaborate in the project.

### Disclaimer 2
If you are going to use OpenBB for your reasearch please cite me: Massimo De Mauri (massimo.demauri@kuleuven.be), KU Leuven, Leuven, Belgium. (and of course cite the authors of the solvers you are going to use)

[1] B. Stellato, G. Banjac, P. Goulart, A. Bemporad and S. Boyd; An Operator Splitting Solver for Quadratic Programs; 2017.
