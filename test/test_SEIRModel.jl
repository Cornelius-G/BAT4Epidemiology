Test.@testset "SEIRModel" begin
    ages = ["A1", "A2"]
    regions = ["R1", "R2", "R3"]

    nages = length(ages)
    nregions = length(regions)
    cshape = (nages=nages, nregions=nregions)

    u0_s = reshape(fill(997, nages*nregions), (nages, nregions))
    u0_e = reshape(fill(2, nages*nregions), (nages, nregions))
    u0_i = reshape(fill(1, nages*nregions), (nages, nregions))
    u0_r = reshape(fill(0., nages*nregions), (nages, nregions))
    u0 = (susceptible=u0_s, exposed=u0_e, infectious=u0_i, recovering=u0_r)

    t_range = 0.:1.:300.

    prior = BAT.NamedTupleDist(
        α = product_distribution(fill(Uniform(0, 0.06), nages*nregions)),
        β = product_distribution(fill(Uniform(1e-5, 40e-5), nages*nregions)),
        γ = product_distribution(fill(Uniform(0, 0.06), nages*nregions)),
    )

    model = SEIRModel(u0, t_range, BAT.varshape(prior), cshape)
    @test isa(model, SEIRModel)


    truth_params = (
        α = fill(0.03, nages*nregions), 
        β = fill(20e-5, nages*nregions), 
        γ = fill(0.03, nages*nregions) 
    )

    sol = solve_ode(model, truth_params)
    @test isa(sol, BAT4Epidemiology.SciMLBase.ODESolution)

    fwm = forward_model(model, truth_params)
    @test isa(fwm, BAT.ValueShapes.NamedTupleDist)

    simulated_data = rand(fwm)
    @test isa(simulated_data, NamedTuple)
    @test simulated_data.t == 0.:1.0:300.
    @test size(simulated_data.new_positive) == (nages, nregions, length(t_range))
    
    loglikelihood = likelihood(model, simulated_data)
    @test isa(loglikelihood, Function)

    posterior = PosteriorDensity(loglikelihood, prior)
    @test isa(posterior, PosteriorDensity)

    vsp = BAT.varshape(prior)
    x = inverse(vsp)(truth_params)
    v_dual = vsp(BAT.ForwardDiff.Dual.(x, zero(x)))

    @inferred BAT.logdensityof(posterior, truth_params)
    @inferred BAT.logdensityof(posterior, v_dual)
    

    alg = MaxDensityLBFGS()
    r_mode = bat_findmode(posterior, alg).result
    @test isa(r_mode, NamedTuple)

    v = truth_params
    f = BAT.logdensityof(posterior)
    gradlogp = BAT.valgradof(f)
    @test isapprox(gradlogp(v)[1], f(v), rtol=0.01)
    
    df_v = BAT.ForwardDiff.gradient(BAT.unshaped(f), (BAT.unshaped(v, BAT.varshape(posterior))))
    gdf_v = BAT.unshaped(gradlogp(v)[2], BAT.varshape(posterior))
    @test isapprox(df_v, gdf_v, rtol=0.01) 
end