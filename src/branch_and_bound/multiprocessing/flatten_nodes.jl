# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-22T16:19:51+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: flatten_nodes.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-11T19:44:41+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}



# this function returns the size of a flat representation of a node
function flat_size(numVars::Int,numDscVars::Int,numCnss::Int)::Int
    return 4 + 4 +        # header + (average fractionality + objective + pseudo-objective + reliable)
           2*numDscVars + # branching bounds
           2*numVars +    # primal + bound_dual
           numCnss        # cns_dual
end


function flat_size(node::BBnode)::Int

    numVars = length(node.primal)
    numDscVars = length(node.branchLoBs)
    numCnss = length(node.cnsDual)

    return flat_size(numVars,numDscVars,numCnss)
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
    numDscVars = length(node.branchLoBs)
    numCnss = length(node.cnsDual)

    @assert length(destinationArray) >= flat_size(numVars,numDscVars,numCnss) + offset

    # header
    destinationArray[offset+1] = 1. # type of node
    destinationArray[offset+2] = numDscVars
    destinationArray[offset+3] = numVars
    destinationArray[offset+4] = numCnss
    offset += 4

    # numeric values
    destinationArray[offset+1] = node.avgAbsFrac
    destinationArray[offset+2] = node.objective
    destinationArray[offset+3] = node.pseudoObjective
    destinationArray[offset+4] = node.reliable
    offset += 4

    # bounds
    @. destinationArray[offset+1:offset+numDscVars] = node.branchLoBs
    offset+=numDscVars
    @. destinationArray[offset+1:offset+numDscVars] = node.branchUpBs
    offset+=numDscVars

    # primal
    @. destinationArray[offset+1:offset+numVars] = node.primal
    offset += numVars

    # dual
    @. destinationArray[offset+1:offset+numVars] = node.bndDual
    offset += numVars
    @. destinationArray[offset+1:offset+numCnss] = node.cnsDual
    offset += numCnss

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
        numDscVars = Int(flatRepresentation[offset+2])
        numVars    = Int(flatRepresentation[offset+3])
        numCnss    = Int(flatRepresentation[offset+4])
        offset += 4

        # numeric values
        avgAbsFrac      = flatRepresentation[offset+1]
        objective       = flatRepresentation[offset+2]
        pseudoObjective = flatRepresentation[offset+3]
        reliable        = Bool(flatRepresentation[offset+4])
        offset += 4

        # bounds
        branchLoBs = Array{Float64,1}(undef,numDscVars)
        @. branchLoBs = flatRepresentation[offset+1:offset+numDscVars]
        offset += numDscVars
        branchUpBs = Array{Float64,1}(undef,numDscVars)
        @. branchUpBs = flatRepresentation[offset+1:offset+numDscVars]
        offset += numDscVars

        # primal
        primal = Array{Float64,1}(undef,numVars)
        @. primal = flatRepresentation[offset+1:offset+numVars]
        offset += numVars

        # dual
        bndDual = Array{Float64,1}(undef,numVars)
        @. bndDual = flatRepresentation[offset+1:offset+numVars]
        offset += numVars
        cnsDual = Array{Float64,1}(undef,numCnss)
        @. cnsDual = flatRepresentation[offset+1:offset+numCnss]
        offset += numCnss

        return BBnode(branchLoBs,branchUpBs,primal,bndDual,cnsDual,avgAbsFrac,objective,pseudoObjective,reliable)
    elseif flatRepresentation[offset+1] == -1.0
        return KillerNode(flatRepresentation[offset+2])
    else
        println(flatRepresentation)
        @error "type of node unknown"
        return
    end
end
