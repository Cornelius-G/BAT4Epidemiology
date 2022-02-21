export SEIRModel
export solve_ode, forward_model


struct SEIRModel <: CompartmentalModel
    u0::NamedTuple
    t_range::AbstractVector{<:Real}
    param_shape::AbstractValueShape
    compartments::NamedTuple
end

compartments(model::SEIRModel) = model.compartments
param_shape(model::SEIRModel) = model.param_shape
state_shape(model::SEIRModel) = BAT.valshape(model.u0)
compartment_shape(model::SEIRModel) = values(compartments(model))
nages(model::SEIRModel) = model.compartments.nages
nregions(model::SEIRModel) = model.compartments.nregions
nsteps(model::SEIRModel) = length(model.t_range)
allshape(model::SEIRModel) = (compartment_shape(model)..., nsteps(model))


function (model::SEIRModel)( # functor to make instances of SEIRModel callable
    du::AbstractVector{<:Real}, 
    u::AbstractVector{<:Real}, 
    p::Any, 
    t::Real
)
    eqn_system(model, du, u, p, t)
end


function eqn_system(model::SEIRModel, du, u , p ,t)
    vs_u = state_shape(model)
    vs_v = param_shape(model)
 
    α = reshape(getindex(p, vs_v.α), compartment_shape(model))
    β = reshape(getindex(p, vs_v.β), compartment_shape(model))
    γ = reshape(getindex(p, vs_v.γ), compartment_shape(model))

    susceptible = view(u, vs_u.susceptible)
    exposed = view(u, vs_u.exposed)
    infectious = view(u, vs_u.infectious)
    recovering = view(u, vs_u.recovering)

    d_susceptible = view(du, vs_u.susceptible)
    d_exposed = view(du, vs_u.exposed)
    d_infectious = view(du, vs_u.infectious)
    d_recovering = view(du, vs_u.recovering)


    for a in 1:nages(model), r in 1:nregions(model)
        d_susceptible[a,r] = -β[a,r]*susceptible[a,r]*infectious[a,r]
        d_exposed[a,r] = β[a,r]*susceptible[a,r]*infectious[a, r] - α[a,r]*exposed[a,r]
        d_infectious[a,r] = α[a,r]*exposed[a,r] - γ[a,r]*infectious[a,r]
        d_recovering[a,r] = γ[a,r]*infectious[a,r]
    end
end


function solve_ode(model::SEIRModel, v; solver=Tsit5(), saveat=model.t_range, dt=1)
    t_min = minimum(model.t_range)
    t_max = maximum(model.t_range)
    v_array = reduce(vcat, v) # use ValueShapes.unshaped ?

    u0 = BAT.unshaped(model.u0, state_shape(model))

    prob = ODEProblem{true}(model, u0, (t_min, t_max), v_array)
    sol = solve(prob, solver, saveat = saveat)

    return sol
end


function generate_observables(model::SEIRModel, v::NamedTuple, U)
    vs_v = param_shape(model)
    α = reshape(v.α, compartment_shape(model))
    #α = reshape(getindex(p, vs_v.α), compartment_shape(model))

    S = flatview(U.susceptible)
    E = flatview(U.exposed)
    I = flatview(U.infectious)
    R = flatview(U.recovering)

    new_positive = Array{eltype(S)}(undef, allshape(model))

    for a in 1:nages(model), r in 1:nregions(model), t in 1:nsteps(model)
       new_positive[a,r,t] = α[a,r] * E[a,r, t]
    end

    return (new_positive = new_positive, )
end


