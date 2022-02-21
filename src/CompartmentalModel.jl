abstract type CompartmentalModel end

export forward_model
export likelihood


shouldbepos(x::Real, tol::Real) = -tol <= x < 0 ? max(x, zero(x)) + eps(typeof(x)) : x


function forward_model(
    model::CompartmentalModel,
    v::NamedTuple, 
)
    sol = solve_ode(model, v)

    result = VectorOfSimilarArrays(sol.u)
    u = state_shape(model).(result)
    #@showme u
    observables_expected = generate_observables(model, v, u)

   return observable_distribution(model, observables_expected)
end


function observable_distribution(
    model::CompartmentalModel,
    obs_expected::NamedTuple#={...}=#,
)
    new_pos = flatview(obs_expected.new_positive)
    new_pos_vec = vec(new_pos)

    dist_new_pos = product_distribution(Poisson.(new_pos_vec))

    return NamedTupleDist(
        t = model.t_range,
        new_positive = reshape(dist_new_pos, allshape(model))
    )
end


function likelihood(model::CompartmentalModel, observed_data::NamedTuple)
    v -> begin
        ll = logpdf(forward_model(model, v), observed_data)
        fixed_ll = isinf(ll) && ll < 0 ? typeof(ll)(-1e10) : ll
        
        return (logval = fixed_ll,)
    end
end
