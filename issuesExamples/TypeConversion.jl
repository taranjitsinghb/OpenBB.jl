# @Author: Wim Van Roy <WimVanRoy>
# @Date:   Fr Sep 27 12:04:25 CEST 2019
# @Email:  wim.vr@hotmail.com
# @Filename: TypeConversion.jl
# @Last modified by:   WimVanRoy
# @Last modified time: Fr Sep 27 12:05:01 CEST 2019
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# Dummy file just to show error type conversions

using SparseArrays
using LinearAlgebra
using OpenBB

obj = OpenBB.LinearObjective(L=zeros(4))
cnss = OpenBB.LinearConstraintSet(A=zeros(4,4), loBs=zeros(4), upBs=zeros(4))
vars = OpenBB.VariableSet(loBs=zeros(4), upBs=zeros(4), vals=zeros(4))

problem = OpenBB.Problem(objFun=obj, cnsSet=cnss, varSet=vars)

workspace = OpenBB.setup(problem, OpenBB.BBsettings(verbose=true), OpenBB.OSQPsettings())
