using Revise, Parameters, Plots
using BifurcationKit

recordFromSolution(x,p) = (u1 = x[1], u2 = x[2])

function pp2!(dz, z, p, t = 0)
    @unpack p1, p2, p3, p4 = p
    u1, u2 = z
    dz[1] = p2 * u1 * (1 - u1) - u1 * u2 - p1 * (1 - exp(-p3 * u1))
    dz[2] = -u2 + p4 *u1 * u2
    dz
end

par_pp2 = (p1 = 1., p2 = 3., p3 = 5., p4 = 3.)

z0 = zeros(2)

prob = BifurcationProblem(pp2!, z0, par_pp2,
    # specify the continuation parameters
    (@lens _.p1), record_from_solution = recordFromSolution)

opts_br = ContinuationPar(p_min = 0.1, p_max = 1.0, dsmax = 0.01,
    #options to detect bifurcations
    detect_bifurcation = 3, n_inversion = 8, max_bisection_steps = 25,
    # number of eigenvalues
    nev = 2,
    # maximum number of contination stepts
    max_steps = 1000,)

diagram = bifurcationdiagram(prob, PALC(),
    3,
    (args...) -> setproperties(opts_br; ds = -0.001, dsmax = 0.01, n_inversion = 8, detect_bifurcation = 3);
    # δp = -0.01
    verbosity = 0, plot = false)

scene = plot(diagram; code = (), title="$(size(diagram)) branches", legend = false)


brH = get_branch(diagram, (2, 1)).γ

optn_po = NewtonPar(tol = 1e-8, max_iterations = 25)

opts_po_cont = ContinuationPar(dsmax = 0.1, ds = -0.001, dsmin = 1e-4, newton_options = (@set optn_po.tol = 1e-8), tol_stability = 1e-2, detect_bifurcation = 1)

Mt = 101 # number of time sections
	br_po = continuation(
	brH, 1, opts_po_cont,
	PeriodicOrbitTrapProblem(M = Mt;
	    # specific linear solver for ODEs
    	jacobian = :Dense);
	record_from_solution = (x, p) -> (xtt=reshape(x[1:end-1],2,Mt);
		return (max = maximum(xtt[1,:]),
			min = minimum(xtt[1,:]),
			period = x[end])),
	finalise_solution = (z, tau, step, contResult; prob = nothing, kwargs...) -> begin
		# limit the period
		z.u[end] < 100
		end,
	normC = norminf)
plot(diagram); plot!(br_po, label = "Periodic orbits", legend = :bottomright)

