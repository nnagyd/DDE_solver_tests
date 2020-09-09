using DifferentialEquations, Plots, DelimitedFiles

function logistic1(dx,x,h,p,t)
    τ = p
    hist = h(p,t-τ)
    dx[1] = -hist[1]
end

#parameters
τ = 1
delays = [τ]
p = τ
tStart = 0.0
tEnd = 10.0
tspan = (tStart,tEnd)
x0 = [1.0]
tol = 13
tolerance = 10.0^(-1*tol)

#initial function
h(p, t) = 1.0

#creating DDE problem
prob = DDEProblem(logistic1,x0,h,tspan,p; constant_lags = delays)
alg = MethodOfSteps(Tsit5())

#solution
sol = solve(prob,alg, dt = 0.01, abstol = tolerance,reltol=tolerance)
plot(sol)

#save points
io = open("tol_"*string(tol)*".txt","w")
writedlm(io,[sol.t collect(Iterators.flatten(sol.u))])
close(io)
