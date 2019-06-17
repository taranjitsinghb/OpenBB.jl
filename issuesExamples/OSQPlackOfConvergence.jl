# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-07T15:30:35+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: OSQPlackOfConvergence.jl
# @Last modified by:   massimo
# @Last modified time: 2019-05-07T15:45:51+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

using MAT
using OpenBB


# open the problem info file
file = matopen(Base.source_dir()*"/problem_matrices.mat")
Q_ = read(file, "P")
L_ = convert(Array{Float64,1},read(file, "q")[:])
A_ = read(file, "A")
cnsLoBs_ = read(file, "lbg")[:]
cnsUpBs_ = read(file, "ubg")[:]
varLoBs_ = read(file, "lbv")[:]
varUpBs_ = read(file, "ubv")[:]
dscIndices_ = convert(Array{Int64, 1},read(file, "discrete_idx")[:]) # Discrete indices (0 or 1)
sos1Groups_ = repeat([-1],length(dscIndices_))

# problem info
numStates = 18
numControls = size(varLoBs_)[1] - 2*numStates
numConstraints = size(A_)[1]
initialState = read(file, "x0")[:]
close(file);

# set parameters
varLoBs_[19] = varUpBs_[19] = 0

# remove info on the first state as it is fixed
varLoBs_ = varLoBs_[numStates+1:end]
varUpBs_ = varUpBs_[numStates+1:end]

# construct the full time horizon problem
horizonLen = 10
# variable values and bounds
varLoBs = vcat(initialState,repeat(varLoBs_,horizonLen))
varUpBs = vcat(initialState,repeat(varUpBs_,horizonLen))
varVals = zeros(size(varLoBs))
# variable values and bounds
varLoBs = vcat(initialState,repeat(varLoBs_,horizonLen))
varUpBs = vcat(initialState,repeat(varUpBs_,horizonLen))
varVals = zeros(size(varLoBs))

# first set of discrete variables
dscIndices = copy(dscIndices_)
sos1Groups = copy(sos1Groups_)
# eliminate the redundant discrete indices
dscToKeep = findall(@. dscIndices_ > numStates)
dscIndices_ = dscIndices_[dscToKeep]
sos1Groups_ = sos1Groups_[dscToKeep]
# discrete variables for the full horizon
for s in 2:horizonLen
    append!(dscIndices,@. dscIndices_ + (s-1)*numStates + (s-1)*numControls)
    append!(sos1Groups,@. sos1Groups_ + (s-1)*(sos1Groups_>-1))
end
vars = OpenBB.VariableSet(loBs=varLoBs,upBs=varUpBs,val=varVals,dscIndices=dscIndices)


# objective
L = L_[1:numStates]
for s in 1:horizonLen
    append!(L, L_[numStates+1:end])
end
objf = OpenBB.LinearObjective(L=L)

# constraints
cnsLoBs = repeat(cnsLoBs_,horizonLen)
cnsUpBs = repeat(cnsUpBs_,horizonLen)
A = sparse(A_)
for s in 2:horizonLen
    global A = vcat(hcat(A,spzeros(numConstraints*(s-1), numControls + numStates)),
                    hcat(spzeros(numConstraints,(s-1) * (numControls + numStates)), A_))
end
cnss = OpenBB.LinearConstraintSet(A=A,loBs=cnsLoBs,upBs=cnsUpBs)

problem = OpenBB.Problem(objFun=objf,cnsSet=cnss,varSet=vars)


if OpenBB.withGUROBI()
    println("solving with gurobi-QP")
    workspace1 = OpenBB.setup(problem,
                              OpenBB.BBsettings(verbose=false,dynamicMode=true,statusInfoPeriod=1),
                              OpenBB.GUROBIsettings())

    OpenBB.solve!(workspace1)
    println(" - solved!")
end

println("solving with OSQP")
workspace2 = OpenBB.setup(problem,
                         OpenBB.BBsettings(verbose=false,dynamicMode=true,statusInfoPeriod=1),
                         OpenBB.OSQPsettings())

OpenBB.solve!(workspace2)
println(" - solved!")
