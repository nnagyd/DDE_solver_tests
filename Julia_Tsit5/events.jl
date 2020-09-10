using DifferentialEquations, Plots, DelimitedFiles, CPUTime, Statistics

const nrOfRuns = 2
const N = 128
const nrOfP = N
const pList = collect(LinRange(0.1,1.5,nrOfP))
const nrOfτ = N
const τList = collect(LinRange(0.1,2.0,nrOfP))
tol = 8
tolerance = 10.0^(-1*tol)
bounces = zeros(nrOfP*nrOfτ,3)
nocheck(dt,u,p,t) = false

function events(dx,x,h,p,t)
    τ1 = p[1]
    hist1 = h(p,t-τ1)
    dx[1] = -9.81
    dx[2] = x[1]-p[2] * hist1[1]
end

#initial function
function h(p,t)
    return [0,20]
end

#event when x[2] == 0
function condition(x,t,integrator)
  x[2]
end

#change sign of x[1]
function affect!(integrator)
  integrator.u[1] = -0.95*integrator.u[1]
end

cb = ContinuousCallback(condition,affect!,save_positions=(false,true),rootfind = true,abstol=1e-6,reltol=1e-6)
alg = MethodOfSteps(Tsit5())

#inital conditions
tstart = 0.0
tend = 60.0
tspan = (tstart,tend)

#solution
times = Vector{Float64}(undef,nrOfRuns)
for runs in 1:nrOfRuns
    tStart = CPUtime_us()
    for i in 1:nrOfτ
        for j in 1:nrOfP
            τ = τList[i]
            p = pList[j]
            #parameters
            delays = [τ]
            pars = [τ,p]
            x0 = [0,20]

            #creating DDE problem
            prob = DDEProblem(events,x0,h,tspan,pars; constant_lags = delays)

            global sol = solve(
                prob,
                alg,
                adaptive = true,
                save_everystep = false,
                save_end = false,
                save_start = false,
                dense = false,
                dt = 0.01,
                abstol = tolerance,
                reltol=tolerance,
                callback = cb,
                unstable_check = nocheck
                )
            bounces[(i-1)*nrOfP + j,1] = p
            bounces[(i-1)*nrOfP + j,2] = τ
            bounces[(i-1)*nrOfP + j,3] = length(sol.t)
        end
    end
    tEnd = CPUtime_us()
    times[runs] = (tEnd-tStart)/(10^6)
    println(times[runs])
end

println("Parameter number: "*string(N*N))
println("Time: "*string(times[1:nrOfRuns]))
println("Avg: "*string(mean(times[2:nrOfRuns])))

#save points
io = open("endvals.txt","w")
writedlm(io,bounces)
close(io)
