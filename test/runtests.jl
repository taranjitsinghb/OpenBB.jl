# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:19:53+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: runtests.jl
# @Last modified by:   massimo
# @Last modified time: 2019-08-26T17:25:33+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

println("OpenBB tests:")
include("./test_problem_definitions.jl")
include("./test_preprocessing.jl")
include("./test_QP_subsolvers.jl")
include("./test_flat_interface.jl")

if OpenBB.withMPCaddon()
    include(Base.source_dir()*"../../../MPCforOpenBB/test/run_tests.jl")
end
