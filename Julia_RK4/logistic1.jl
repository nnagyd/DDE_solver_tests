using DifferentialEquations, Plots, DelimitedFiles, LoopVectorization, CPUTime, Statistics

const unroll = 256
const steps = 10000
const nrOfRuns = 5
const nrOfParameters = 8192
const parameterList = collect(LinRange(0.0,4.0,nrOfParameters))
nocheck(dt,u,p,t) = false


function logistic1(dx,x,h,p,t)
    τ = p[1]
    hist = h(p,t-τ)
    @avx for i in 1:unroll
        dx[i] = x[i]*(p[i+1]-hist[i])
    end
end

#delays
τ = 1.0
delays = τ * ones(unroll)

#parameters
p = Vector{Float64}(undef,unroll + 1)
p[1] = τ
for i in 2:unroll+1
    p[i] = 1.0
end

#initial condition
x0 = 0.5*ones(unroll)

#time
tspan = (0.0,10.0)
dt = 10.0 / steps

#initial function
h(p, t) = ones(unroll) * (1.5 - cos(t))

#solver algorithm
alg = MethodOfSteps(RK4())

#creating DDE problem
prob = DDEProblem(logistic1,x0,h,tspan,p; constant_lags = delays)

#ensemble parameter changhe and output function
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
io = open("julia_endvals.txt","w")
writedlm(io,[parameterList collect(Iterators.flatten(ensSol.u))])
close(io)
