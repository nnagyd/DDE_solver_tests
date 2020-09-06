using DifferentialEquations, Plots, DelimitedFiles, CPUTime, Statistics

const nrOfRuns = 1
const unroll = 128 #unrolls the inner loop, where p changes
const N = 128
const nrOfP = N
const pList = collect(LinRange(0.1,1.5,nrOfP))
const nrOfτ = N
const τList = collect(LinRange(0.1,2.0,nrOfP))
const nrOfSteps =  10000
bounces = zeros(nrOfP*nrOfτ,3)
counter = zeros(unroll)

function events(dx,x,h,p,t)
    τ1 = p[1]
    h(hist,p,t-τ1)
    for i in 1:unroll
        ofs = 2*(i-1)
        dx[ofs+1] = -9.81
        dx[ofs+2] = x[ofs+1]-p[i+1] * hist[ofs+1]
    end
end

#initial function
hist = zeros(2*unroll)
function h(hist,p,t)
    for i in 1:2:unroll*2
        hist[i] = 0
        hist[i+1] = 20
    end
end

#event when x[2] == 0
function condition(out,x,t,integrator)
  for i in 1:unroll
      out[i] = x[2*i]
  end
end

#change sign of x[1]
function affect!(integrator, idx)
    for i in 1:unroll
        if idx == i
            ofs = 2*i - 1
            integrator.u[ofs] = -0.95*integrator.u[ofs]
            global counter[idx] += 1
        end
    end
end

cb = VectorContinuousCallback(condition,affect!,unroll,save_positions=(false,true),rootfind = true,abstol=1e-6,reltol=0)
alg = MethodOfSteps(RK4())

#inital conditions
tstart = 0.0
tend = 60.0
tspan = (tstart,tend)
dt = tend / nrOfSteps

#solution
times = Vector{Float64}(undef,nrOfRuns)
for runs in 1:nrOfRuns
    tStart = CPUtime_us()
    for i in 1:nrOfτ
        for j in 1:unroll:nrOfP
            τ = τList[i]

            #parameters
            delays = [τ]
            pars = zeros(unroll+1)
            pars[1] = τ
            pars[2:unroll+1] = pList[j:j+unroll-1]

            x0 = zeros(unroll*2)
            for i in 1:2:unroll*2
                x0[i] = 0
                x0[i+1] = 20
            end

            #creating DDE problem
            prob = DDEProblem(events,x0,h,tspan,pars; constant_lags = delays)
            global counter = zeros(unroll)

            global sol = solve(
                prob,
                alg,
                adaptive = false,
                save_everystep = false,
                save_end = false,
                save_start = false,
                dense = false,
                dt = dt,
                callback = cb
                )
            for k in 1:unroll
                bounces[(i-1)*nrOfP + j + k - 1,1] = pList[j+k-1]
                bounces[(i-1)*nrOfP + j + k - 1,2] = τ
                bounces[(i-1)*nrOfP + j + k - 1,3] = counter[k]
            end

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
