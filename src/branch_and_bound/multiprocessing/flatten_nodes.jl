# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-22T16:19:51+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: flatten_nodes.jl
# @Last modified by:   massimo
# @Last modified time: 2019-08-27T20:10:43+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}



# this function returns the size of a flat representation of a node
function flat_size(numVars::Int,numCnss::Int)::Int
    return 3 +            # header
           4 +            # average fractionality + objective + pseudo-objective + reliable
           4*numVars +    # variable bounds + primal + bound_dual
           3*numCnss      # constraints bounds + constraints_dual
end


function flat_size(node::BBnode)::Int
    return flat_size(length(node.primal),length(node.cnsDual))
end

function flat_size(node::NullBBnode)::Int
    return 1
end

function flat_size(node::KillerNode)::Int
    return 2
end

# fill the destination array with the flat representation of the given node
function flatten_in!(node::BBnode,destinationArray::T;offset::Int=0)::Int where T <: AbstractArray

    numVars = length(node.primal)
    numCnss = length(node.cnsDual)

    @assert length(destinationArray) >= flat_size(numVars,numCnss) + offset

    # header
    destinationArray[offset+1] = 1. # type of node
    destinationArray[offset+2] = numVars
    destinationArray[offset+3] = numCnss
    offset += 3

    # numeric values
    destinationArray[offset+1] = node.avgAbsFrac
    destinationArray[offset+2] = node.objective
    destinationArray[offset+3] = node.pseudoObjective
    destinationArray[offset+4] = node.reliable
    offset += 4

    # bounds
    @. destinationArray[offset+1:offset+numVars] = node.varLoBs; offset+=numVars
    @. destinationArray[offset+1:offset+numVars] = node.varUpBs; offset+=numVars
    @. destinationArray[offset+1:offset+numCnss] = node.cnsLoBs; offset+=numCnss
    @. destinationArray[offset+1:offset+numCnss] = node.cnsUpBs; offset+=numCnss

    # primal
    @. destinationArray[offset+1:offset+numVars] = node.primal; offset += numVars

    # dual
    @. destinationArray[offset+1:offset+numVars] = node.bndDual; offset += numVars
    @. destinationArray[offset+1:offset+numCnss] = node.cnsDual; offset += numCnss

    return offset

end

function flatten_in!(node::NullBBnode,destinationArray::T;offset::Int=0)::Int where T <: AbstractArray
    @assert length(destinationArray) >= 1 + offset
    destinationArray[offset+1] = 0. # type of node
    offset += 1
    return offset
end

function flatten_in!(node::KillerNode,destinationArray::T;offset::Int=0)::Int where T <: AbstractArray
    @assert length(destinationArray) >= 2 + offset
    destinationArray[offset+1] = -1. # type of node
    destinationArray[offset+2] = node.count
    offset += 2
    return offset
end


# this function construct and returns a flat array representation of a BBnode
function flatten(node::BBnode;offset::Int=0)::Array{Float64,1}

    # initialize target array
    flatRepresentation = zeros(flat_size(node)+offset)
    flatten_in!(node,flatRepresentation,offset=offset)

    return flatRepresentation
end



# reconstruct a node from its flat representation
function rebuild_node(flatRepresentation::T1;offset::Int=0)::AbstractBBnode where T1 <: AbstractArray


    # read type
    if  flatRepresentation[offset+1] == 0.0
        return NullBBnode()

    elseif flatRepresentation[offset+1] == 1.0

        # rest of header
        numVars = Int(flatRepresentation[offset+2])
        numCnss = Int(flatRepresentation[offset+3])
        offset += 3

        # numeric values
        avgAbsFrac      = flatRepresentation[offset+1]
        objective       = flatRepresentation[offset+2]
        pseudoObjective = flatRepresentation[offset+3]
        reliable        = Bool(flatRepresentation[offset+4])
        offset += 4

        # bounds
        varLoBs = flatRepresentation[offset+1:offset+numVars]; offset += numVars
        varUpBs = flatRepresentation[offset+1:offset+numVars]; offset += numVars
        cnsLoBs = flatRepresentation[offset+1:offset+cnsVars]; offset += cnsVars
        cnsUpBs = flatRepresentation[offset+1:offset+cnsVars]; offset += cnsVars

        # primal
        primal = flatRepresentation[offset+1:offset+numVars]; offset += numVars

        # dual
        bndDual = flatRepresentation[offset+1:offset+numVars]; offset += numVars
        cnsDual = flatRepresentation[offset+1:offset+numCnss]; offset += numCnss

        return BBnode(varLoBs,varUpBs,cnsLoBs,cnsUpBs,primal,bndDual,cnsDual,avgAbsFrac,objective,pseudoObjective,reliable)

    elseif flatRepresentation[offset+1] == -1.0
        return KillerNode(flatRepresentation[offset+2])
    else
        println(flatRepresentation)
        @error "type of node unknown"
        return
    end
end
