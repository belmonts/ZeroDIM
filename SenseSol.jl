"""
A Julia code meant to perform Global Sensisitivity Analysis for Zero Dimensional Biofilm

Author:
    - Alexandra Mazanko (@belmontz)

"""
#############################################
# Last Update: September 20, 2023
#############################################
#Include all the necessary packages
using LinearAlgebra
using DifferentialEquations
using Plots
using DiffEqSensitivity
using QuasiMonteCarlo
using GlobalSensitivity

function InitialConditions()

    u0 = [0.5; 0.5; 0.5; 0.5]

    return u0
end

function GlobalSenseSolu(tspan::Tuple, u0::Vector)

    function reactionparamter_definition()
        q = 1.0
        μ = 6.0
        δ = 6.0
        α = 6.0
        Y = 2.0
        β = 1.0
        ks = 4.0
        Sn = 25.0
        rpvc = [q,μ,δ, α, Y, β, ks, Sn]
        return rpvc
    end

    function RHSfun(du, u, p, t)
        parms  = zeros(8)
        parms[1:8] = rpd

        parms[parmsToVary] = p

        rpvc = parms[1:8]

        q = rpvc[1]
        μ = rpvc[2]
        δ = rpvc[3]
        α = rpvc[4]
        Y = rpvc[5]
        β = rpvc[6]
        ks = rpvc[7]
        Sn = rpvc[8]
      
        du[1] = -q*u[1] + μ*(u[3]/(ks + u[3]))*u[1] + δ*u[2] - α*u[1]
        du[2] = μ*(u[3]/(ks + u[3]))*u[2] - δ*u[2] + α*u[1]
        du[3] = q*(Sn - u[3]) - (μ/Y)*(u[3]/(ks + u[3]))*(u[1] + u[2])
        du[4] = -q*u[4] + β*(u[3]/(ks + u[3]))*(u[1] + u[2])
    end

    u0 = InitialConditions()

    rpvc = reactionparamter_definition()
    parm = zeros(8)
    parm[1:8] = rpvc*1.5
    #parm[1:8] = [rpvc[1:2];rpvc[3:4]*10;rpvc[5:8]]

    global parmsToVary = 1:8

    global rpd = reactionparamter_definition()

    prob = ODEProblem(RHSfun, u0, tspan, parm)
    t = collect(range(0, stop=10, length=100))
    alg = Rosenbrock23()

    N = 100000
    lb = parm
    ub = parmsToVary
    sampler = SobolSample()
    A,B = QuasiMonteCarlo.generate_design_matrices(N,lb,ub,sampler)

    f1 = function (p)
        prob1 = remake(prob;p=p)
        sol  = solve(prob1, alg; saveat=t)
        [sol[end][3], sol[end][4]]#we make our designated region the final solution here.
    end

     sobol_result = GlobalSensitivity.gsa(f1, Sobol(), A,B)
  
    return sobol_result

end
