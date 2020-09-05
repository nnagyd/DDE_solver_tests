using DifferentialEquations, Plots, DelimitedFiles,LoopVectorization, CPUTime, Statistics

include("C:/Users/nnagy/Documents/Egyetem/HDS/DDE/Julia/create_mesh.jl")

const unroll = 1024
const steps = 10000
const nrOfRuns = 3
const nrOfParameters = 262144
const parameterList = collect(LinRange(0.0,2.0,nrOfParameters))
nocheck(dt,u,p,t) = false

function logistic2(dx,x,h,p,t)
    τ1 = p[1]
    τ2 = p[2]
    hist1 = h(p,t-τ1)
    hist2 = h(p,t-τ2)
    @avx for i in 1:unroll
        dx[i] = x[i]*hist2[i]*(p[i+2]-hist1[i])
    end
end


#delays
delays = Vector{Float64}(undef,unroll*2)
τ1 = 1.0
τ2 = 2.0
for i in 1:unroll
    ofs = 2*(i-1)
    delays[ofs + 1] = τ1
    delays[ofs + 2] = τ2
end

#parameters
p = Vector{Float64}(undef,unroll + 2)
p[1] = τ1
p[2] = τ2
for i in 3:unroll+2
    p[i] = 1.0
end

#initial conditions
x0 = ones(unroll)*0.5

#time
tstart = 0.0
tend = 10.0
tspan = (tstart,tend)
dt = tend / steps

#mesh
disc = [-1.5,-√2,-1.101,-0.5]
mesh = finalMesh(disc,[τ1,τ2],tStart,tEnd)

#initial function
function h(p,t)
    if t < -1.5
        return cos(4 * pi * t) * ones(unroll)
    elseif t < -sqrt(2)
        return t * t * ones(unroll)
    elseif t < -1.101
        return exp(t) * ones(unroll)
    elseif t < -0.5
        return 0 * ones(unroll)
    else
        return (t + 0.5) * ones(unroll)
    end
end

#solver algorithm
alg = MethodOfSteps(RK4())

#creating DDE problem
prob = DDEProblem(logistic2,x0,h,tspan,p; constant_lags = delays)

#ensemble parameter changhe and output function
function newPar(prob,i,repeat)
    ofs = unroll*(i-1)
    for j in 1:unroll
        prob.p[j+2] = parameterList[ofs+j]
    end
    if i % 100 == 0
        println(i)
    end
    prob
end

outputSave(sol,i) = (sol[end],false)

#creating ensemble problem
ensProb = EnsembleProblem(prob,prob_func = newPar,output_func = outputSave)

#solution
times = Vector{Float64}(undef,nrOfRuns)
for runs in 1:nrOfRuns
    tStart = CPUtime_us()
    global ensSol = solve(
        ensProb,
        alg,
        EnsembleSerial(),
        adaptive = false,
        dt = dt,
        trajectories = nrOfParameters ÷ unroll,
        save_everystep = false,
        dense = false,
        unstable_check = nocheck,
        d_discontinuities = mesh
        )
    tEnd = CPUtime_us()
    times[runs] = (tEnd-tStart)/(10^6)
    println(times[runs])
end

println("Parameter number: "*string(nrOfParameters))
println("Time: "*string(times[1:nrOfRuns]))
println("Avg: "*string(mean(times[2:nrOfRuns])))

#save points
io = open("C:/Users/nnagy/Documents/Egyetem/HDS/DDE/results/julia_logistic2/julia_endvals.txt","w")
writedlm(io,[parameterList collect(Iterators.flatten(ensSol.u))])
close(io)
