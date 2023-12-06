#Author: belmontz, sept 5th, 2023

#Call all the packages we will be using for our sensitivity analysis
using Plots
using DifferentialEquations
using LinearAlgebra
using GlobalSensitivity
pyplot()

#Define our function containing all the reaction parameters of our problem
function reactionparameter_def()
    q = 1.2
    μ = 6.0
    δ = 1.0
    α = 0.01
    Y = 2.0
    β = 1.0
    ks = 4.0
    Sn = 25.0
    rpvc = [q,μ,δ, α, Y, β, ks, Sn]
    return rpvc
end

function monod(c, ks)

  if c <= 0.0
    return m = 0.0
  else
    return m = c/(c + ks)
  end

  return m
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

    du[1] = μ*(u[2] /(ks + u[2]))*u[1] - q*u[1]
    du[2] = q*(Sn - u[2]) - y*(μ*(u[2] /(ks + u[2])))*u[1]
    du[3] = β*(μ*(u[2] /(ks + u[2])))*u[1] - q*u[3]
  end

#Globally define our timescale, initial conditions, and reaction parameters for CSTR
#u0 = [0.0; 25.0; 0.0]

rpvc = reactionparameter_def()

#Define our given problem
#prob = ODEProblem(simpleCSTR, u0, tspan, rpvc)

#sol = solve(prob)

function zerodim(du, u, p, t)
  q = p[1]
  μ = p[2]
  δ = p[3]
  α = p[4]
  Y = p[5]
  β = p[6]
  ks = p[7]
  Sn = p[8]

  du[1] = -q*u[1] + μ*monod(u[3], ks)*u[1] - δ*u[2] + α*u[1]
  du[2] = μ*(monod(u[3], ks))*u[2] + δ*u[2] - α*u[1]
  du[3] = q*(Sn - u[3]) - (μ/Y)*(monod(u[3], ks))*(u[1] + u[2])
  du[4] = -q*u[4] + β*(monod(u[3], ks))*(u[1] + u[2])
end



function zerodim2D(du, u, p, t)
  q = p[1]
  μ = p[2]
  δ = p[3]
  α = p[4]
  Y = p[5]
  β = p[6]
  ks = p[7]
  Sn = p[8]

  du[1] = -q*u[1] + μ*monod(u[3], ks)*u[1] - δ*u[2]^2 + α*u[1]^2
  du[2] = μ*(monod(u[3], ks))*u[2] + δ*u[2]^2 - α*u[1]^2
  du[3] = q*(Sn - u[3]) - (μ/Y)*(monod(u[3], ks))*(u[1] + u[2])
  du[4] = -q*u[4] + β*(monod(u[3], ks))*(u[1] + u[2])
end


u0 = [1.0; 0.0; 5.0; 0.0]

tspan = (0.0, 10.0)

prob = ODEProblem(zerodim, u0, tspan, rpvc)

prob2D = ODEProblem(zerodim2D, u0, tspan, rpvc)

sol = solve(prob)

sol2D = solve(prob2D)


dim1plot = plot(sol)

dim2plot = plot(sol2D)

plot(dim1plot, dim2plot)



