export moving_average
export average
export small_sol

moving_average(vs,n) = [sum(@view vs[i:(i+n-1)])/n for i in 1:(length(vs)-(n-1))]


function average(data::NamedTuple, avg_span::Real) # for NamedTuple with shape (t =, new_positive =)
    new_positive_avg = moving_average(collect(data.new_positive), avg_span)
    t = 0.:1.:(length(new_positive_avg)-1)

    return (t = t, new_positive = new_positive_avg)
end


function small_sol(sol, dim)
    u = sol.u 
    n = Int(length(u[1])/4)
    l = size(u, 1)

    S = [u[i][dim] for i in 1:l]
    E = [u[i][n+dim] for i in 1:l]
    I = [u[i][2n+dim] for i in 1:l]
    R = [u[i][3n+dim] for i in 1:l]
    return (t=sol.t, S=S, E=E, I=I, R=R)
end