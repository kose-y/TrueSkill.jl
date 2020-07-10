mutable struct TrueSkillEnv
    mu::Float64
    sigma::Float64
    beta::Float64
    tau::Float64
    draw_probability::Float64
end

GLOBAL_ENV = TrueSkillEnv(25.0, 25.0/3, 25.0/6, 25.0/300, 0.10)
DELTA = 0.0001

function calc_draw_probability(draw_margin, size; env=GLOBAL_ENV)
    return 2 * cdf(draw_margin / (sqrt(size) * env.beta)) - 1.
end

function calc_draw_margin(draw_probability, size; env=GLOBAL_ENV)
    return ppf((draw_probability + 1) / 2.) * sqrt(size) * env.beta
end

function _team_sizes(rating_groups)
    szs = map(x -> length(x), rating_groups)
    cumsum(szs)
end

function _floating_point_error()
    error("Floating point error. Result out of range.")
end

struct Rating <: GaussianItem
    pi::Float64
    tau::Float64
end

function Rating(; mu=GLOBAL_ENV.mu, sigma=GLOBAL_ENV.sigma)
    pi = sigma^(-2)
    tau = pi * mu
    Rating(pi, tau)
end

function Rating(g::Gaussian)
    Rating(g.pi, g.tau)
end

@inline function mu(r::Rating)
    if r.pi != 0
        return r.tau / r.pi
    else
        return 0
    end
end

@inline function sigma(r::Rating)
    if r.pi != 0
        return(sqrt(1/r.pi))
    else
        return Inf
    end
end

@inline function isless(a::Rating, b::Rating)
    isless(mu(a), mu(b))
end

function show(io::IO, r::Rating)
    print(io, "Rating(mu=$(mu(r)), sigma=$(sigma(r)))")
end

function create_rating(;mu=GLOBAL_ENV.mu, sigma=GLOBAL_ENV.sigma)
    Rating(;mu=mu, sigma=sigma)
end

function v_win(diff::Real, draw_margin::Real)
    x = diff - draw_margin
    denom = cdf(x)
    if denom != 0.0
        return pdf(x) / denom
    else
        return -x
    end
end

function v_draw(diff::Real, draw_margin::Real)
    abs_diff = abs(diff)
    a, b = draw_margin - abs_diff, -draw_margin - abs_diff
    denom = cdf(a) - cdf(b)
    numer = pdf(b) - pdf(a)
    c = 0.0
    if denom != 0.0
        c = numer / denom
    else
        c = a
    end
    d = 0
    if diff < 0
        d = -1
    else
        d = 1
    end
    return c * d
end

function w_win(diff::Real, draw_margin::Real)
    x = diff - draw_margin
    v = v_win(diff, draw_margin)
    w = v * (v + x)
    if 0 < w < 1
        return w
    else
        _floating_point_error()
    end
end

function w_draw(diff::Real, draw_margin::Real)
    abs_diff = abs(diff)
    a, b = draw_margin - abs_diff, -draw_margin - abs_diff
    denom = cdf(a) - cdf(b)
    if denom == 0
        _floating_point_error()
    end
    v = v_draw(abs_diff, draw_margin)
    return (v ^ 2) + (a * pdf(a) - b * pdf(b)) / denom
end

function validate_rating_groups(rating_groups::Vector{Vector{Rating}})
    if length(rating_groups) < 2
        error("need multiple groups")
    end
    if any(map(length, rating_groups) .== 0)
        error("no empty groups allowed")
    end
    return rating_groups
end

function factor_graph_builders(rating_groups::Vector{Vector{Rating}}, ranks, weights; env=GLOBAL_ENV)
    # This is the main builder!!
    flatten_ratings = reduce(vcat, rating_groups)
    flatten_weights = reduce(vcat, weights)
    size = length(flatten_ratings)
    group_size = length(rating_groups)

    rating_vars = [Variable() for x in 1:size]
    perf_vars = [Variable() for x in 1:size]
    team_perf_vars = [Variable() for x in 1:group_size]
    team_diff_vars = [Variable() for x in 1:(group_size - 1)]
    team_sizes = _team_sizes(rating_groups)

    # layer builders
    function build_rating_layer()
        r = Vector{PriorFactor}()
        for (rating_var, rating) in zip(rating_vars, flatten_ratings)
            push!(r, PriorFactor(rating_var, rating, env.tau))
        end
        return r
    end

    function build_perf_layer()
        r = Vector{LikelihoodFactor}()
        for (rating_var, perf_var) in zip(rating_vars, perf_vars)
            push!(r, LikelihoodFactor(rating_var, perf_var, env.beta ^ 2))
        end
        return r
    end

    function build_team_perf_layer()
        r = Vector{SumFactor}()
        for (team, team_perf_var) in enumerate(team_perf_vars)
            if team > 1
                first = team_sizes[team - 1] + 1
            else
                first = 1
            end
            last = team_sizes[team]
            child_perf_vars = perf_vars[first:last]
            coeffs = flatten_weights[first:last]
            push!(r, SumFactor(team_perf_var, child_perf_vars, coeffs))
        end
        return r
    end

    function build_team_diff_layer()
        r = Vector{SumFactor}()
        for (team, team_diff_var) in enumerate(team_diff_vars)
            push!(r, SumFactor(team_diff_var, 
                team_perf_vars[team:team+1], [+1, -1]))
        end
        return r
    end

    function build_trunc_layer() 
        r = Vector{TruncateFactor}()
        for (x, team_diff_var) in enumerate(team_diff_vars)
            draw_probability = env.draw_probability
            size = sum(map(length, rating_groups[x:x+1]))
            draw_margin = calc_draw_margin(draw_probability, size; env=env)
            if ranks[x] == ranks[x + 1] # in case of tie
                v_func, w_func = v_draw, w_draw
            else
                v_func, w_func = v_win, w_win
            end
            push!(r, TruncateFactor(team_diff_var, v_func, w_func, draw_margin))
        end
        return r
    end
    return (build_rating_layer, build_perf_layer, build_team_perf_layer,
        build_team_diff_layer, build_trunc_layer)
end

function run_schedule(build_rating_layer, build_perf_layer, 
    build_team_perf_layer, build_team_diff_layer, build_trunc_layer; min_delta=DELTA)
    # This is the main runner!!
    if min_delta <= 0
        error("min_delta must be greater than 0")
    end
    layers = Vector{Vector{Node}}()

    rating_layer = build_rating_layer()
    perf_layer = build_perf_layer()
    team_perf_layer = build_team_perf_layer()

    push!(layers, rating_layer)
    push!(layers, perf_layer)
    push!(layers, team_perf_layer)

    for f in reduce(vcat, layers)
        down(f)
    end

    team_diff_layer = build_team_diff_layer()
    trunc_layer = build_trunc_layer()

    push!(layers, team_diff_layer)
    push!(layers, trunc_layer)

    team_diff_len = length(team_diff_layer)

    for i in 1:10
        if team_diff_len == 1
            # two teams
            down(team_diff_layer[1])
            delta = up(trunc_layer[1])
        else
            # three or more teams
            delta = 0.0
            for j in 1:(team_diff_len - 1)
                down(team_diff_layer[j])
                delta = max(delta, up(trunc_layer[j]))
                up(team_diff_layer[j]; index=2) # up to right Variable
            end
            for j in team_diff_len:-1:2
                down(team_diff_layer[j])
                delta = max(delta, up(trunc_layer[j]))
                up(team_diff_layer[j]; index=1) # up to left Variable
            end
        end
        if delta <= min_delta
            break
        end
    end

    # up both ends
    up(team_diff_layer[1]; index=1)
    up(team_diff_layer[team_diff_len]; index=2)

    # up the remainder of the black arrows
    for f in team_perf_layer
        for i in 1:(length(f.vars) - 1) # number of members.
            up(f; index=i)
        end
    end
    for f in perf_layer
        up(f)
    end
    return layers
end

function rate(rating_groups::Vector{Vector{Rating}}; ranks=1:length(rating_groups), weights=[ones(length(x)) for x in rating_groups], min_delta=DELTA, env=GLOBAL_ENV)
    validate_rating_groups(rating_groups)
    group_size = length(rating_groups)
    if length(ranks) != group_size
        error("wrong ranks")
    end
    
    # sort rating groups by ranks
    by_rank = x -> x[2][2]
    sorting = sort(collect(enumerate(zip(rating_groups, ranks, weights))), by=by_rank)
    sorted_rating_groups = Vector{Vector{Rating}}()
    sorted_ranks = Vector{Float64}()
    sorted_weights = Vector{Vector{Float64}}()
    for (x, (g, r, w)) in sorting
        push!(sorted_rating_groups, g)
        push!(sorted_ranks, r)
        push!(sorted_weights, [max(min_delta, w_) for w_ in w])
    end
    b1, b2, b3, b4, b5 = factor_graph_builders(sorted_rating_groups, sorted_ranks, sorted_weights; env=env)
    layers = run_schedule(b1, b2, b3, b4, b5; min_delta=min_delta)

    rating_layer, team_sizes = layers[1], _team_sizes(sorted_rating_groups)
    transformed_groups = Vector{Vector{Rating}}()
    for (firstm1, last) in zip(vcat([0], team_sizes[1:end-1]), team_sizes)
        group = Vector{Rating}()
        for f in rating_layer[(firstm1+1):last]
            push!(group, Rating(; mu=mu(var(f)), sigma=sigma(var(f))))
        end
        push!(transformed_groups, group)
    end

    by_hint = x -> x[1]
    unsorting = sort(collect(zip((x for (x, _) in sorting), transformed_groups)), by=by_hint)
    return [g for (x, g) in unsorting]
end

function expose(rating::Rating; env::TrueSkillEnv=GLOBAL_ENV)
    k = env.mu / env.sigma
    return mu(rating) - k * sigma(rating)
end

function make_as_global(env::TrueSkillEnv)
    GLOBAL_ENV = env
end

function rate(rating1::Rating, rating2::Rating; drawn=false, min_delta=DELTA)
    if drawn
        ranks = [0, 0]
    else
        ranks = [0 ,1]
    end
    teams = rate([[rating1], [rating2]]; ranks=ranks,  min_delta = min_delta)
    return teams[1][1], teams[2][1]
end

function rate(ratings::Vector{Rating}; ranks=1:length(ratings), min_delta=DELTA) 
    teams = rate(map(x -> [x], ratings); ranks=ranks, min_delta=min_delta)
    return map(x -> x[1], teams)
end
