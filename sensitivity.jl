#Author: belmontz, sept 5th, 2023

#Call all the packages we will be using for our sensitivity analysis
using Plots
using DifferentialEquations
using LinearAlgebra
using GlobalSensitivity
pyplot()

#Define our function containing all the reaction parameters of our problem
function reactionparameter_def()
    q = 0.01
    μ = 6.0
    δ = 0.01
    α = 0.01
    Y = 2.0
    β = 1.0
    ks = 4.0
    Sn = 1000.0
    rpvc = [q,μ,δ, α, Y, β, ks, Sn]
    return rpvc
end

#Define our single species CSTR model
function simpleCSTR(du, u, p, t)
    q = p[1]
    μ = p[2]
    δ = p[3]
    α = p[4]
    Y = p[5]
    β = p[6]
    ks = p[7]
    Sn = p[8]

    y = 1/Y

    du[1] = (μ*(u[2] /(ks + u[2]))*u[1] - q*u[1]
    du[2] = q*(Sn - u[2]) - y*(μ*(u[2] /(ks + u[2])))*u[1]
    du[3] = β*(μ*(u[2] /(ks + u[2])))*u[1] - q*u[3]
  end

#Globally define our timescale, initial conditions, and reaction parameters for CSTR
tspan = (0.0, 25.0)
u0 = [0.0; 25.0; 0.0]
rpvc = reactionparameter_def()

#Define our given problem
prob = ODEProblem(simpleCSTR, u0, tspan, rpvc)

sol = solve(prob)

animate(sol)
