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
u0 = [0.0; 25.0; 0.0]
rpvc = reactionparameter_def()

#Define our given problem
prob = ODEProblem(simpleCSTR, u0, tspan, rpvc)

sol = solve(prob)

function zerodim(du, u, p, t)
  q = p[1]
  μ = p[2]
  δ = p[3]
  α = p[4]
  Y = p[5]
  β = p[6]
  ks = p[7]
  Sn = p[8]

  du[1] = -q*u[1] + μ*(u[3]/(ks + u[3]))*u[1] + δ*u[2] - α*u[1]
  du[2] = μ*(u[3]/(ks + u[3]))*u[2] - δ*u[2] + α*u[1]
  du[3] = q*(Sn - u[3]) - (μ/Y)*(u[3]/(ks + u[3]))*(u[1] + u[2])
  du[4] = -q*u[4] + β*(u[3]/(ks + u[3]))*(u[1] + u[2])
end

u0z = [1.0; 1.0; 0.0; 0.0]

probz = ODEProblem(zerodim, u0z, tspan, rpvc)

solz = solve(probz)

#Local sensitivity problem




animate(sol)
t = collect(range(0, stop=10, length=100))
alg = Rosenbrock23()

f1 = function (p)
  prob1 = remake(probz;p=p)
  sol  = solve(prob1, alg; saveat=t)
  [sol[end][3]]#we make our designated region the final solution here.
end

#We must convert our rpvc vector into a tuple for sensitivity analysis

low_bound = rpvc*0.95

low_bound = tuple(low_bound...)

up_bound  = rpvc*1.05

up_bound = tuple(up_bound...)




param_range = [low_bound, up_bound]

sobol_result = GlobalSensitivity.gsa(f1,Sobol(), [low_bound, up_bound], samples = 1000)



