# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-22T16:19:51+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: flatten_nodes.jl
# @Last modified by:   massimo
# @Last modified time: 2019-05-31T13:27:59+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}



# this function returns the size of a flat representation of a node
function flat_size(numVars::Int,numDscVars::Int,numCnss::Int)::Int
    return 4 + 3 +        # header + average fractionality + objective value + reliable
           5*numDscVars + # branching bounds + pseudoCosts
           2*numVars +    # primal + bound_dual
           numCnss        # cns_dual
end


function flat_size(node::BBnode)::Int

    numVars = length(node.primal)
    numDscVars = length(node.pseudoCosts)
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
    numDscVars = length(node.pseudoCosts)
    numCnss = length(node.cnsDual)

    @assert length(destinationArray) >= flat_size(numVars,numDscVars,numCnss) + offset

    # header
    destinationArray[offset+1] = 1. # type of node
    destinationArray[offset+2] = numDscVars
    destinationArray[offset+3] = numVars
    destinationArray[offset+4] = numCnss
    offset += 4

    # numeric values
    destinationArray[offset+1] = node.objVal
    destinationArray[offset+2] = node.avgFrac
    destinationArray[offset+3] = node.reliable
    offset += 3

    # bounds
    for (k,pair) in enumerate(node.branchLoBs)
        destinationArray[offset+k] = pair[1]
        destinationArray[offset+numDscVars+k] = pair[2]
    end
    if length(node.branchLoBs) < numDscVars # put a zero to mark the end of the info
        destinationArray[offset+length(node.branchLoBs)+1] = 0.
    end
    offset+=2*numDscVars
    for (k,pair) in enumerate(node.branchUpBs)
        destinationArray[offset+k] = pair[1]
        destinationArray[offset+numDscVars+k] = pair[2]
    end
    if length(node.branchUpBs) < numDscVars # put a zero to mark the end of the info
        destinationArray[offset+length(node.branchUpBs)+1] = 0.
    end
    offset+=2*numDscVars

    # pseudo-costs
    @. destinationArray[offset+1:offset+numDscVars] = node.pseudoCosts
    offset += numDscVars

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
        numVars =    Int(flatRepresentation[offset+3])
        numCnss =    Int(flatRepresentation[offset+4])
        offset += 4

        # numeric values
        objVal = flatRepresentation[offset+1]
        avgFrac = flatRepresentation[offset+2]
        reliable = Bool(flatRepresentation[offset+ 3])
        offset += 3

        # bounds
        branchLoBs = Dict{Int,Float64}()
        for k in 1:numDscVars
            if flatRepresentation[offset+k] != 0.
                branchLoBs[Int(flatRepresentation[offset+k])] = flatRepresentation[offset + numDscVars + k]
            else
                break
            end
        end
        offset += 2*numDscVars

        branchUpBs = Dict{Int,Float64}()
        for k in 1:numDscVars
            if flatRepresentation[offset+k] != 0.
                branchUpBs[Int(flatRepresentation[offset+k])] = flatRepresentation[offset + numDscVars + k]
            else
                break
            end
        end
        offset += 2*numDscVars

        #pseudo costs
        pseudoCosts = Array{Float64,1}(undef,numDscVars)
        @. pseudoCosts = flatRepresentation[offset+1:offset+numDscVars]
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

        return BBnode(branchLoBs,branchUpBs,pseudoCosts,primal,bndDual,cnsDual,avgFrac,objVal,reliable)
    elseif flatRepresentation[offset+1] == -1.0
        return KillerNode(flatRepresentation[offset+2])
    else
        println(flatRepresentation)
        @error "type of node unknown"
        return
    end
end
