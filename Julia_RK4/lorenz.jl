using DifferentialEquations, DelimitedFiles, LoopVectorization, CPUTime, Statistics

const unroll = 512
const steps = 10000
const nrOfRuns = 3
const nrOfParameters = 8192
const parameterList = collect(LinRange(0.0,50.0,nrOfParameters))
nocheck(dt,u,p,t) = false

function lorenz(dx,x,h,p,t)
    τ = p[1]
    h(hist,p,t-τ)
    @avx for i in 1:unroll
        ofs = 3 * (i-1)
        dx[ofs+1] = 10.0*(hist[ofs+2] - x[ofs+1])
        dx[ofs+2] = p[i+1]*x[ofs+1]-x[ofs+2]-x[ofs+3]*x[ofs+1]
        dx[ofs+3] = x[ofs+1]*x[ofs+2] - (8.0/3.0)* x[ofs+3]
    end
end

#delays
τ = 13.0/28.0
delays = τ * ones(unroll)

#parameters
p = Vector{Float64}(undef,unroll + 1)
p[1] = τ
for i in 2:unroll+1
    p[i] = 50.0
end

#initial conditions
x0 = -8*ones(unroll*3)

#time
tstart = 0.0
tend = 10.0
tspan = (tstart,tend)
dt = tend / steps

#initial function
const hist = zeros(unroll*3)
function h(hist,p,t)
    for i = 1:3:unroll*3
        hist[i] = -8
        hist[i+1]=-8+sin(2*pi*t)
        hist[i+2]=-8
    end
end

#solver algorithm
alg = MethodOfSteps(RK4())

#creating DDE problem
prob = DDEProblem(lorenz,x0,h,tspan,p; constant_lags = delays)

#ensemble parameter change and output function
function newPar(prob,i,repeat)
    ofs = unroll*(i-1)
    for j in 1:unroll
        prob.p[j+1] = parameterList[ofs+j]
    end
    if i % 100 == 0
        println(i)
    end
    prob
end

outputSave(sol,i) = (sol[end],false)

#creating ensemble problem
ensProb = EnsembleProblem(prob,prob_func = newPar) #output_func = outputSave

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
        unstable_check = nocheck
        )

    tEnd = CPUtime_us()
    times[runs] = (tEnd-tStart)/(10^6)
    println(times[runs])
    GC.gc()
end

println("Parameter number: "*string(nrOfParameters))
println("Time: "*string(times[1:nrOfRuns]))
println("Avg: "*string(mean(times[2:nrOfRuns])))

#save points
io = open("julia_endvals.txt","w")
endvals = zeros(nrOfParameters)
idx = 1
for sol in ensSol
    for i in 1:3:unroll*3
        endvals[idx] = sol[end][i]
        global idx += 1
    end
end
writedlm(io,[parameterList endvals])
close(io)
