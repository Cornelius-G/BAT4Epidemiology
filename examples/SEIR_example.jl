using Pkg
Pkg.activate(abspath(joinpath(@__DIR__, "..")))
Pkg.instantiate()

using BAT4Epidemiology
using BAT
using Distributions
using Plots
using AxisKeys
using InverseFunctions
using ShowMe

gr(size=(1.1*850, 1.1*600), thickness_scaling = 1.3)
plot = Plots.plot

include("../src/CovidData.jl")

#====== SEIR model with generic truth data ====================#
ages = ["A1", "A2"]
regions = ["R1", "R2", "R3"]
nages = length(ages)
nregions = length(regions)
cshape = (nages=nages, nregions=nregions)

u0_s = reshape(fill(997, nages*nregions), (nages, nregions))
u0_e = reshape(fill(2, nages*nregions), (nages, nregions))
u0_i = reshape(fill(1, nages*nregions), (nages, nregions))
u0_r = reshape(fill(0., nages*nregions), (nages, nregions))

A = KeyedArray(u0_s; ages=ages, regions=regions)
B = KeyedArray(u0_e; ages=ages, regions=regions)
C = KeyedArray(u0_i; ages=ages, regions=regions)
D = KeyedArray(u0_r; ages=ages, regions=regions)
u0 = (susceptible=A, exposed=B, infectious=C, recovering=D)

t_range = 0.:1.:300.

prior = BAT.NamedTupleDist(
    α = product_distribution(fill(Uniform(0, 0.06), nages*nregions)),
    β = product_distribution(fill(Uniform(1e-5, 40e-5), nages*nregions)),
    γ = product_distribution(fill(Uniform(0, 0.06), nages*nregions)),
)

model = SEIRModel(u0, t_range, BAT.varshape(prior), cshape)

truth_params = rand(prior)
vsp = BAT.varshape(prior)
x = inverse(vsp)(truth_params)
v_dual=vsp(BAT.ForwardDiff.Dual.(x, zero(x)))

# solve & plot truth ODEs
sol = solve_ode(model, truth_params)
#sol = solve_ode(model, v_dual)
#plot(sol, label=["Succeptible" "Exposed" "Infected" "Removed"])

# simulate fake data
fwm = forward_model(model, truth_params)
fwm = forward_model(model, v_dual)

simulated_data = rand(fwm)
plot(simulated_data.t, simulated_data.new_positive[1, 3, :], yguide = "new positive", xguide = "days", label = "simulated data")

#testing moving average
#simulated_data_averaged = average(simulated_data, 5) 
#plot!(simulated_data_averaged.t, simulated_data_averaged.new_positive[1, 3, :])

# Bayesian inference on simulated data
loglikelihood = likelihood(model, simulated_data)
posterior = PosteriorDensity(loglikelihood, prior)

BAT.logdensityof(posterior, v_dual)
BAT.logdensityof(posterior, truth_params)

#sampling_algorithm = MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^5, nchains = 4)
sampling_algorithm = MCMCSampling(mcalg = HamiltonianMC())
samples = bat_sample(posterior, sampling_algorithm).result
plot(samples)

# find mode of posterior distribution
alg = MaxDensityLBFGS()
#alg = MaxDensityNelderMead()
r_mode = bat_findmode(posterior, alg).result

# plot best fit together with (simulated) data and truth model
plot(t_range, mean(fwm.new_positive), label="best fit")
plot!(simulated_data.t, simulated_data.new_positive[:], yguide = "new positive", xguide = "days", label = "simulated data")
plot!(simulated_data_averaged.t,simulated_data_averaged.new_positive, label = "simulated averaged")
