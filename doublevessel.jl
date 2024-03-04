#Author: belmontz, Mar 04, 2024
#Call all the packages we will be using for our sensitivity analysis
using DifferentialEquations
using LinearAlgebra
using CairoMakie

#Define our reaction parameters
function reactionparameter_def()
    q = 1.2
    ν1 = 5.9
    ν2 = 5.9
    μ1 = 6.0
    μ2 = 6.1
    D1 = 1.0
    D2 = 2.5
    So = 25.0
    ks = 0.7
    rpvc = [q, ν1, ν2, μ1, μ2, D1, D2, So, ks]
    return rpvc
end

rpvc = reactionparameter_def()

function monod(c, ks)

    if c <= 0.0
      return m = 0.0
    else
      return m = c/(c + ks)
    end
  
    return m
  end

function dualvessel(du, u, p,t)
    q = p[1]
    ν1 = p[2]
    ν2 = p[3]
    μ1 = p[4]
    μ2 = p[5]
    D1 = p[6]
    D2 = p[7]
    So = p[8]
    ks = p[9]
    du[1] = q*(So - u[1]) - ν1*monod(u[1],ks)*u[3] - D1*u[3] + D2*u[4]
    du[2] = -ν2*monod(u[2],ks)*u[4] - D2*u[4] + D1*u[3]
    du[3] = μ1*monod(u[1],ks)*u[3] - q*u[3]  - D1*u[3]  +D2*u[4]
    du[4] = μ2*monod(u[2],ks)*u[4] + D1*u[3] - D2*u[4]

end

u0 = [1.0; 0.0; 5.0; 0.0]

tspan = (0.0, 10.0)

prob = ODEProblem(dualvessel, u0, tspan, rpvc)
sol = solve(prob)
using Plots
Plots.plot(sol)