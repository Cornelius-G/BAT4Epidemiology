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
using Test

gr(size=(1.1*850, 1.1*600), thickness_scaling = 1.3)
plot = Plots.plot

include("../src/CovidData.jl")

#====== SEIR model with generic truth data ====================#
ages = ["A1",]
regions = ["R1",]
nages = length(ages)
nregions = length(regions)
cshape = (nages=nages, nregions=nregions)

u0_s = reshape(fill(9997, nages*nregions), (nages, nregions))
u0_e = reshape(fill(20, nages*nregions), (nages, nregions))
u0_i = reshape(fill(10, nages*nregions), (nages, nregions))
u0_r = reshape(fill(0., nages*nregions), (nages, nregions))

A = KeyedArray(u0_s; ages=ages, regions=regions)
B = KeyedArray(u0_e; ages=ages, regions=regions)
C = KeyedArray(u0_i; ages=ages, regions=regions)
D = KeyedArray(u0_r; ages=ages, regions=regions)
u0 = (susceptible=A, exposed=B, infectious=C, recovering=D)

t_range = 1.:1.:28.

# nchangepoints = 4 # number of changepoints
# lndagger_pd = product_distribution(fill(Normal(4, 1), nages*nregions*nchangepoints))
# lndagger = reshape(lndagger_pd, (nages, nregions, nchangepoints) )

# Δdn_pd = product_distribution(fill(Normal(0, 3.5), nages*nregions*nchangepoints))
# Δdn = reshape(Δdn_pd, (nages, nregions, nchangepoints))

# Δγ_pd = product_distribution(fill(Normal(0, 0.5), nages*nregions*nchangepoints))
# Δγ = reshape(Δγ_pd, (nages, nregions, nchangepoints))

# prior = BAT.NamedTupleDist(
#     α = product_distribution(fill(Uniform(0, 0.6), nages*nregions)),
#     γ = product_distribution(fill(Uniform(0, 0.05), nages*nregions)),
#     R0 = product_distribution(fill(truncated(Normal(0.01, 0.000001), 0, Inf), nages*nregions)),
#     Δγ = Δγ_pd,
#     lndagger = lndagger_pd,
#     Δdn = Δdn_pd,
# )


nchangepoints = 4 # number of changepoints
lndagger_pd = fill(4, nages*nregions*nchangepoints)
lndagger = reshape(lndagger_pd, (nages, nregions, nchangepoints) )

Δdn_pd = fill(0., nages*nregions*nchangepoints)
Δdn = reshape(Δdn_pd, (nages, nregions, nchangepoints))

Δγ_pd = fill(0.1, nages*nregions*nchangepoints)
Δγ = reshape(Δγ_pd, (nages, nregions, nchangepoints))


prior = BAT.NamedTupleDist(
    α = product_distribution(fill(Uniform(0, 0.6), nages*nregions)),
    γ = product_distribution(fill(Uniform(0, 0.05), nages*nregions)),
    R0 = fill(0.01, nages*nregions),
    Δγ = Δγ_pd,
    lndagger = lndagger_pd,
    Δdn = Δdn_pd,
)
BAT.vardof(prior)

model = CPSEIRModel(u0, t_range, BAT.varshape(prior), cshape, nchangepoints)

truth_params = rand(prior)

vsp = BAT.varshape(prior)
x = inverse(vsp)(truth_params)
v_dual=vsp(BAT.ForwardDiff.Dual.(x, zero(x)))


# solve & plot truth ODEs
sol = solve_ode(model, truth_params)
s_sol = small_sol(sol, 1)

plot(s_sol.t, s_sol.S, label="Susceptible")
plot!(s_sol.t, s_sol.E, label="Exposed")
plot!(s_sol.t, s_sol.I, label="Infectious")
plot!(s_sol.t, s_sol.R, label="Recovered")


# simulate fake data
fwm_dual = forward_model(model, v_dual)
fwm = forward_model(model, truth_params)

simulated_data = rand(fwm)
plot(simulated_data.t, simulated_data.new_positive[1, 1, :], yguide = "new positive", xguide = "days", label = "simulated data")

#testing moving average
#simulated_data_averaged = average(simulated_data, 5) 
#plot!(simulated_data_averaged.t, simulated_data_averaged.new_positive[1, 1, :])

# Bayesian inference on simulated data
llikelihood = likelihood(model, simulated_data)
posterior = PosteriorDensity(llikelihood, prior)
@inferred BAT.logdensityof(posterior, v_dual)
@inferred BAT.logdensityof(posterior, truth_params)

sampling_algorithm = MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^5, nchains = 4)
#sampling_algorithm = MCMCSampling(mcalg = HamiltonianMC())
samples = bat_sample(posterior, sampling_algorithm).result
#plot(samples)

# find mode of posterior distribution
alg = MaxDensityLBFGS()
#alg = MaxDensityNelderMead()

r_mode = bat_findmode(posterior, alg).result
truth_params

r_mode_flat = inverse(vsp)(r_mode)
truth_params_flat = inverse(vsp)(truth_params)

plot(1:length(r_mode_flat), (r_mode_flat .- truth_params_flat)./truth_params_flat, legend=false)

# plot best fit together with (simulated) data and truth model
# plot(t_range, mean(fwm.new_positive), label="best fit")
# plot!(simulated_data.t, simulated_data.new_positive[:], yguide = "new positive", xguide = "days", label = "simulated data")
# plot!(simulated_data_averaged.t,simulated_data_averaged.new_positive, label = "simulated averaged")
