export CPSEIRModel
export solve_ode, forward_model


struct CPSEIRModel <: CompartmentalModel
    u0::NamedTuple
    t_range::AbstractVector{<:Real}
    param_shape::AbstractValueShape
    compartments::NamedTuple
    nchangepoints::Int64
end

compartments(model::CPSEIRModel) = model.compartments
param_shape(model::CPSEIRModel) = model.param_shape
state_shape(model::CPSEIRModel) = BAT.valshape(model.u0)
compartment_shape(model::CPSEIRModel) = values(compartments(model))
nages(model::CPSEIRModel) = model.compartments.nages
nregions(model::CPSEIRModel) = model.compartments.nregions
nsteps(model::CPSEIRModel) = length(model.t_range)
nchangepoints(model::CPSEIRModel) = model.nchangepoints
allshape(model::CPSEIRModel) = (compartment_shape(model)..., nsteps(model))


function (model::CPSEIRModel)( # functor to make instances of CPSEIRModel callable
    du::AbstractVector{<:Real}, 
    u::AbstractVector{<:Real}, 
    p::Any, 
    t::Real
)
    eqn_system(model, du, u, p, t)
end


function Rt(t, R0, Δγ, lndagger, Δdn)
    n = length(Δγ)

    γ = zeros(eltype(Δγ), n)
    for i in 1:n
        γ[i] = γn(t, i-1, Δγ[i], Δdn[i], lndagger[i])
    end

    return R0 * exp(sum(γ))
end


function γn(t, n, Δγn, Δdn, lndagger)
    ln = log(1 + exp(lndagger))
    dn = 7 * n + Δdn

    denom = 1 + exp((-4/ln)*(t-dn))
    return Δγn / denom
end


function eqn_system(model::CPSEIRModel, du, u , p ,t)
    vs_u = state_shape(model)
    vs_v = param_shape(model)

   # @showme vs_u
   # @showme vs_v
   # @showme compartment_shape(model)
 
    α = reshape(getindex(p, vs_v.α), compartment_shape(model))
    γ = reshape(getindex(p, vs_v.γ), compartment_shape(model))
    R0 = reshape(getindex(p, vs_v.R0), compartment_shape(model))
    
    Δγ = reshape(getindex(p, vs_v.Δγ), (compartment_shape(model)..., 4))

    lndagger = reshape(getindex(p, vs_v.lndagger), (compartment_shape(model)..., 4))
    Δdn = reshape(getindex(p, vs_v.Δdn), (compartment_shape(model)..., 4))

    R = ones(eltype(R0), compartment_shape(model))

    for a in 1:nages(model), r in 1:nregions(model)
        R[a, r] = Rt(t, R0[a, r], Δγ[a, r, :], lndagger[a, r, :], Δdn[a, r, :]) 
    end
    #@showstop R[1, 1]

    #--- debugging -------------------
    # t = 1:300
    # a=1; r=1
    # R2 = Rt.(t, Ref(R0[a, r]), Ref(Δγ[a, r, :]), Ref(lndagger[a, r, :]), Ref(Δdn[a, r, :])) 
    # @showme R2
    # β2 = γ[a, r] .* R2
    # @showme β2
    #---------------------------------

    β = γ .* R
   
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


function solve_ode(model::CPSEIRModel, v; solver=Tsit5(), saveat=model.t_range, dt=1)
    t_min = minimum(model.t_range)
    t_max = maximum(model.t_range)
  
    #v_array = reduce(vcat, v) # use ValueShapes.unshaped ?
    vsp = param_shape(model)
    v_array=inverse(vsp)(v)

    u0 = BAT.unshaped(model.u0, state_shape(model))
    

    prob = ODEProblem{true}(model, u0, (t_min, t_max), v_array)
    
    sol = solve(prob, solver, saveat = saveat)

    return sol
end


function generate_observables(model::CPSEIRModel, v::NamedTuple, U)
    vs_v = param_shape(model)
    α = reshape(v.α, compartment_shape(model))
    #α = reshape(getindex(v, vs_v.α), compartment_shape(model))

    S = flatview(U.susceptible)
    E = flatview(U.exposed)
    I = flatview(U.infectious)
    R = flatview(U.recovering)

    new_positive = Array{eltype(S)}(undef, allshape(model))

    for a in axes(new_positive, 1), r in axes(new_positive, 2), t in axes(new_positive, 3)
       new_positive[a,r,t] = α[a,r] * E[a,r, t]
    end

    return (new_positive = new_positive, )
end


