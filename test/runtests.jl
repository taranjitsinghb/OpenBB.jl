# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:19:53+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: runtests.jl
# @Last modified by:   massimo
# @Last modified time: 2019-08-27T19:59:11+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}
using OpenBB

println("OpenBB tests:")
include("./test_problem_definition_fundamentals.jl")
include("./test_problem_definition_update.jl")
include("./test_preprocessing.jl")
include("./test_flat_interface.jl")
include("./test_QP_subsolvers.jl")

if OpenBB.withMPCaddon()
    include(Base.source_dir()*"../../../MPCforOpenBB/test/run_tests.jl")
end
