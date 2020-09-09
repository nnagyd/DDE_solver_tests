using DifferentialEquations, Plots, DelimitedFiles

function analitic2(dx,x,h,p,t)
    τ = p
    hist = h(p,t-τ)
    dx[1] = -hist[1]
end

#parameters
τ = 3
delays = [τ]
p = τ
tStart = 0.0
tEnd = 25.0
tspan = (tStart,tEnd)
x0 = [0.0]
tol = 13
tolerance = 10.0^(-1*tol)

#initial function
function h(p, t)
    if t <= -2
        return 0
    elseif t <= -1
        return 1
    else
        return 0
    end
end

#creating DDE problem
prob = DDEProblem(analitic2,x0,h,tspan,p; constant_lags = delays)
alg = MethodOfSteps(Tsit5())

#solution
sol = solve(prob,alg, dt = 0.01, abstol = tolerance,reltol=tolerance, dtmax = 1.0, dtmin=1e-10)
plot(sol)

#save points
io = open("tol_"*string(tol)*".txt","w")
writedlm(io,[sol.t collect(Iterators.flatten(sol.u))])
close(io)
