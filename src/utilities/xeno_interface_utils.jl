# @Author: Massimo De Mauri <massimo>
# @Date:   2019-06-20T14:09:12+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: xeno_interface_utils.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-20T14:18:41+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


function eval_string(str::String)
    return eval(Base.Meta.parse(str))
end
