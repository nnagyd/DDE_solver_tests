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
steps = 960*4
dt = tEnd / steps

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
alg = MethodOfSteps(RK4())

#solution
sol = solve(prob,alg, adaptive = false, dt = dt)
plot(sol)

#save points
io = open("step_"*string(steps)*".txt","w")
writedlm(io,[sol.t collect(Iterators.flatten(sol.u))])
close(io)
