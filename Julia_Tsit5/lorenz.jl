using DifferentialEquations, DelimitedFiles, CPUTime, Statistics

#number of parameters
const nrOfRuns = 3
const nrOfParameters = 1024
const parameterList = collect(LinRange(0.0,50.0,nrOfParameters))
nocheck(dt,u,p,t) = false

function lorenz(dx,x,h,p,t)
    @inbounds begin
        hist1 = h(p,t-(13.0/28.0))
        dx[1] = 10.0*(hist1[2] - x[1])
        dx[2] = p[1]*x[1]-x[2]-x[3]*x[1]
        dx[3] = x[1]*x[2] - (8.0/3.0)* x[3]
    end
end

#parameters
τ1 = 13.0/28.0
delays = [τ1]
p = [50.0]
tspan = (0.0,10.0)
x0 = [-8,-8,-8]
tol = 8
tolerance = 10.0^(-1*tol)

#initial function
function h(p,t)
    return [-8 -8 + sin(2*pi*t) -8]
end

#solver algorithm
alg = MethodOfSteps(Tsit5())

#creating DDE problem
prob = DDEProblem(lorenz,x0,h,tspan,p; constant_lags = delays)

#ensemble parameter change and output function
function newPar(prob,i,repeat)
    @inbounds begin
        prob.p[1] = parameterList[i]
        if i % 500 == 0
            println(i)
        end
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
        #EnsembleThreads(),
        EnsembleSerial(),
        adaptive = true,
        dt = 0.01,
        trajectories = nrOfParameters,
        save_everystep = false,
        dense = false,
        abstol = tolerance,
        reltol = tolerance,
        unstable_check = nocheck
        )

    tEnd = CPUtime_us()
    times[runs] = (tEnd-tStart)/(10^6)
    println(times[runs])
end

println("Parameter number: "*string(nrOfParameters))
println("Time: "*string(times[1:nrOfRuns]))
println("Avg: "*string(mean(times[2:nrOfRuns])))

#save points
io = open("endvals.txt","w")
solu = permutedims(reshape(hcat(ensSol.u...), (length(ensSol.u[1]), length(ensSol.u))))
writedlm(io,[parameterList collect(Iterators.flatten(solu[:,1]))])
close(io)
