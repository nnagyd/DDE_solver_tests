
function meshDouble(disc, delay, tol, depth)
    discLen = length(disc)
    delayLen = length(delay)
    mesh = Vector{Float64}(undef,discLen*delayLen*2)

    id = 1
    for ξ in disc
        for τ in delay
            mesh[id] = ξ + τ - tol
            id += 1
            mesh[id] = ξ + τ + tol
            id += 1
        end
    end
    return mesh
end

function meshSimple(disc, delay, tol, depth)
    discLen = length(disc)
    delayLen = length(delay)
    mesh = Vector{Float64}(undef,discLen + discLen*delayLen)
    mesh[1:discLen] = disc

    id = discLen+1
    for ξ in disc
        for τ in delay
            mesh[id] = ξ + τ
            id += 1
        end
    end

    if depth == 1
        return mesh
    else
        meshSimple(mesh,delay,tol,depth-1)
    end
end

function finalMesh(disc, delay, tmin, tmax, tol = 1e-10, depth = 4)
    doubleMesh = meshDouble(disc,delay,tol,1)
    removeFromSimple = meshSimple(disc,delay,tol,1)
    simpleMesh = meshSimple(disc,delay,tol,depth)
    simplePoints = setdiff(simpleMesh,removeFromSimple)
    mesh =vcat(simplePoints, doubleMesh)
    sort!(mesh)
    filter!(t-> t > tmin && t < tmax,mesh)
end
