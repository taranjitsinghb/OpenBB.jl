# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:19:45+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: nodes_priority_functions.jl
# @Last modified by:   massimo
# @Last modified time: 2019-03-18T00:31:21+01:00
# @License: apache 2.0
# @Copyright: {{copyright}}


# priority functions for nodes
function lower_objective(l::BBnode,r::BBnode,status::BBstatus)::Bool
    if l.objVal <= r.objVal
        return true
    else
        return false
    end
end

function higher_objective(l::BBnode,r::BBnode,status::BBstatus)::Bool
    if l.objVal >= r.objVal
        return true
    else
        return false
    end
end

function lower_fractionality(l::BBnode,r::BBnode,status::BBstatus)::Bool
    if l.avgFrac <= r.avgFrac
        return true
    else
        return false
    end
end

function higher_fractionality(l::BBnode,r::BBnode,status::BBstatus)::Bool
    if l.avgFrac >= r.avgFrac
        return true
    else
        return false
    end
end

function lower_mixed(l::BBnode,r::BBnode,status::BBstatus)::Bool
    if status.objLoB > -Inf
        if l.objVal + l.avgFrac*(l.objVal-status.objLoB) < r.objVal + r.avgFrac*(r.objVal - status.objLoB)
            return true
        else
            return false
        end

    elseif l.avgFrac <= r.avgFrac
        return true
    else
        return false
    end
end
