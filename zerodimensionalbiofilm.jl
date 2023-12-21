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
    δ = 0.01
    α = 2.0
    Y = 2.0
    β = 1.0
    ks = 4.0
    Sn = 25.0
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

function monod(c, ks)

  if c < 0.0
    return m = 0.0
  else
    m  = c/(ks + c)
  end
end


  #monod function that takes the concentration and half saturation and returns the computed value

  

function zerodim(du, u, p, t)
  q = p[1]
  μ = p[2]
  δ = p[3]
  α = p[4]
  Y = p[5]
  β = p[6]
  ks = p[7]
  Sn = p[8]

  du[1] = -q*u[1] + μ*(u[3]/(ks + u[3]))*u[1] - δ*u[2] + α*u[1]
  du[2] = μ*(u[3]/(ks + u[3]))*u[2] + δ*u[2] - α*u[1]
  du[3] = q*(Sn - u[3]) - (μ/Y)*(u[3]/(ks + u[3]))*(u[1] + u[2])
  du[4] = -q*u[4] + β*(u[3]/(ks + u[3]))*(u[1] + u[2])
end

u0 = [1.0; 0.0; 5.0; 0.0]

prob = ODEProblem(zerodim, u0, tspan, rpvc)

sol = solve(prob)



function zerodimS2(du, u, p, t)
  q = p[1]
  μ = p[2]
  δ = p[3]
  α = p[4]
  Y = p[5]
  β = p[6]
  ks = p[7]
  Sn = p[8]

  du[1] = -q*u[1] + μ*(u[3]/(ks + u[3]))*u[1] + δ*u[2]^2 - α*u[1]^2
  du[2] = μ*(u[3]/(ks + u[3]))*u[2] - δ*u[2]^2 + α*u[1]^2
  du[3] = q*(Sn - u[3]) - (μ/Y)*(u[3]/(ks + u[3]))*(u[1] + u[2])
  du[4] = -q*u[4] + β*(u[3]/(ks + u[3]))*(u[1] + u[2])
end

prob2 = ODEProblem(zerodimS2, u0, tspan, rpvc)

sol2 = solve(prob2)

dim1plot = plot(sol)
dim2plot = plot(sol2)

plot(dim1plot, dim2plot, layout = (1,2))

