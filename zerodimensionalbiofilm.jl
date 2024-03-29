#Author: belmontz, Sept 5th, 2023

#Call all the packages we will be using for our sensitivity analysis
using DifferentialEquations
using LinearAlgebra
using CairoMakie

#Define our function containing all the reaction parameters of our problem
function reactionparameter_def()
    q = 1.2
    μ = 6.0
    δ = 1.1
    α = 1.1
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

  du[1] = -q*u[1] + μ*monod(u[3], ks)*u[1] + δ*u[2] - α*u[1]
  du[2] = μ*(monod(u[3], ks))*u[2] - δ*u[2] + α*u[1]
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

  du[1] = -q*u[1] + μ*monod(u[3], ks)*u[1] + δ*u[2]^2 - α*u[1]^2
  du[2] = μ*(monod(u[3], ks))*u[2] - δ*u[2]^2 + α*u[1]^2
  du[3] = q*(Sn - u[3]) - (μ/Y)*(monod(u[3], ks))*(u[1] + u[2])
  du[4] = -q*u[4] + β*(monod(u[3], ks))*(u[1] + u[2])
end


u0 = [1.0; 0.0; 5.0; 0.0]

tspan = (0.0, 10.0)


function plotter(pv, Numloops, scaler)
  Numloops = Numloops
  psl = Array{Any}(nothing, Numloops)
  rpvc = reactionparameter_def()
  for i in 1:Numloops
    rpvc[pv] = rpvc[pv]*scaler
    lprob = ODEProblem(zerodim, u0, tspan, rpvc)
    lsol = solve(lprob)
    psl[i] = plot(lsol)
  end
  return psl
end

#function that will take two paramters and scale them at the same scalar.
function multiplotter(p1,p2, Numloops, scaler)
  Numloops = Numloops
  psl = Array{Any}(nothing, Numloops)
  rpvc = reactionparameter_def()
  for i in 1:Numloops
    rpvc[p1] = rpvc[p1]*scaler
    rpvc[p2] = rpvc[p2]*scaler
    lprob = ODEProblem(zerodim, u0, tspan, rpvc)
    lsol = solve(lprob)
    psl[i] = plot(lsol)
  end
  return psl
end

function multiplot_stored(p1,p2,Numloops,scaler)
  stl = Array{Any}(nothing,Numloops)
  rpvc = reactionparameter_def()
  for i in 1:Numloops
    rpvc[p1] = rpvc[p1]*scaler
    rpvc[p2] = rpvc[p2]*scaler
    lprob = ODEProblem(zerodim, u0, tspan, rpvc)
    lsol = solve(lprob)
    stl[i] = [rpvc[p1],rpvc[p2],lsol[3][end]]
  end

  #Define our axes
  mx = Array{Any}(nothing,Numloops)
  my = Array{Any}(nothing,Numloops)
  mz = Array{Any}(nothing,Numloops)
  for i in 1:Numloops
    mx[i] = stl[i][1]
    my[i] = stl[i][2]
    mz[i] = stl[i][3]
  end

  mx = Float64.(mx)
  my = Float64.(my)
  mz = Float64.(mz)
  return [mx,my,mz]
end

function multiplot_storedV2(p1,p2,p1_range,p2_range,N_i, N_j)
  stl = Array{Float64,2}(undef,N_i, N_j)
  rpvc = reactionparameter_def()
  mx = Array{Float64}(undef,N_i)
  my = Array{Float64}(undef,N_j)

  for i in 1:N_i
    for j in 1:N_j
      rpvc[p1] = p1_range[2] + (i-1)*((p1_range[1] - p1_range[2])/(N_i - 1))
      rpvc[p2] = p2_range[2]+ (j-1)*((p2_range[1] - p2_range[2])/(N_j - 1))
      lprob = ODEProblem(zerodim, u0, tspan, rpvc)
      lsol = solve(lprob)
      mx[i] = rpvc[p1]
      my[j] = rpvc[p2]
      stl[i,j] = lsol[3][end]
    end
  end

  #Define our axes
  return [mx,my,stl]
end

function multiplot_storedV3(p1,p2,p1_range,p2_range,N_i, N_j,FUN,s)
  stl = Array{Float64,2}(undef,N_i, N_j)
  rpvc = reactionparameter_def()
  mx = Array{Float64}(undef,N_i)
  my = Array{Float64}(undef,N_j)

  for i in 1:N_i
    for j in 1:N_j
      rpvc[p1] = p1_range[2] + (i-1)*((p1_range[1] - p1_range[2])/(N_i - 1))
      rpvc[p2] = p2_range[2]+ (j-1)*((p2_range[1] - p2_range[2])/(N_j - 1))
      lprob = ODEProblem(FUN, u0, tspan, rpvc)
      lsol = solve(lprob)
      mx[i] = rpvc[p1]
      my[j] = rpvc[p2]
      stl[i,j] = lsol[s][end]
    end
  end

  #Define our axes
  return [mx,my,stl]
end

function heatmapper(heatdata)
  fig = Figure()
  hx = heatdata[1]
  hy = heatdata[2]
  hz = heatdata[3]

  ax = Axis(fig[1,2], xlabel = "Detachment Parameter δ", ylabel = "Attachment Parameter α")
  hm = heatmap!(ax,hx,hy,hz)
  Colorbar(fig[1,2],hm)
  
  return fig
end

function heatmapper_flip(heatdata)
  f = Figure()
  hx = heatdata[1]
  hy = heatdata[2]
  hz = heatdata[3]

  heatmap = CairoMakie.heatmap(hx,hy,hz)

  return heatmap
end


function worker_ant(p1, p2,Numloops,scaler)
  work = multiplot_stored(p1,p2,Numloops,scaler)
  heatmap = heatmapper(work)
  return heatmap
end


function twoplotter(p1, p2, NumLoops, s1, s2)
  twopsl = Array{Any}(nothing, NumLoops)
  rpvc = reactionparameter_def()

  if Merge == "True"
    for i in 1:NumLoops
      rpvc[p1] = rpvc[p1]*s1
      rpvc[p2] = rpvc[p2]*s1
      lprob = ODEProblem(zerodim, u0, tspan, rpvc)
      lsol = solve(lprob)
      twopsl[i] = plot(lsol)
    end
  else
    for i in 1:NumLoops
      rpvc[p1] = rpvc[p1]*s1
      rpvc[p2] = rpvc[p2]*s2
      lprob = ODEProblem(zerodim, u0, tspan, rpvc)
      lsol = solve(lprob)
      twopsl[i] = plot(lsol) 
    end
  end
  return twopsl
end

