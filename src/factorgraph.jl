import Base: setindex!, getindex

abstract type Node end


struct Variable <: Node
    val::Gaussian
    messages::Dict{Node, Gaussian}
end

function Variable(; mu::Float64=0.0, sigma::Float64=0.0)
    if sigma == 0
        return Variable(Gaussian(0.0, 0.0), Dict())
    end
    pi = sigma^(-2)
    tau = pi * mu
    Variable(Gaussian(pi, tau), Dict())
end

@inline (*)(a::Variable, b::Variable) = Gaussian(a.val.pi + b.val.pi, a.val.tau + b.val.tau)
@inline (*)(a::Variable, b::Gaussian) = Gaussian(a.val.pi + b.pi, a.val.tau + b.tau)
@inline (*)(a::Gaussian, b::Variable) = Gaussian(a.pi + b.val.pi, a.tau + b.val.tau)
@inline (/)(a::Variable, b::Variable) = Gaussian(a.val.pi - b.val.pi, a.val.tau - b.val.tau)
@inline (/)(a::Variable, b::Gaussian) = Gaussian(a.val.pi - b.pi, a.val.tau - b.tau)
@inline (/)(a::Gaussian, b::Variable) = Gaussian(a.pi - b.val.pi, a.tau - b.val.tau)

@inline function mu(r::Variable)
    if r.val.pi != 0
        return r.val.tau / r.val.pi
    else
        return 0
    end
end

@inline function sigma(r::Variable)
    if r.val.pi != 0
        return(sqrt(1/r.val.pi))
    else
        return Inf
    end
end

function setindex!(v::Variable, message, factor::Node)
    v.messages[factor] = message
end
function getindex(v::Variable, factor::Node)
    v.messages[factor]
end

function delta(v::Variable, w::Gaussian)
    pi_delta = abs(v.val.pi - w.pi)
    if pi_delta == Inf
        return 0.
    end
    return max(abs(v.val.tau - w.tau), sqrt(pi_delta))
end

function delta(v::Variable, w::Variable)
    delta(v, w.val)
end

function set(v::Variable, w::Gaussian)
    _delta = delta(v, w)
    v.val.pi, v.val.tau = w.pi, w.tau
    return _delta
end

function set(v::Variable, w::Variable)
    set(v, w.val)
end

function update_message(v::Variable, factor::Node; pi=0, tau=0, message::Union{Nothing, Gaussian}=nothing)
    if message === nothing
        message = Gaussian(pi, tau)
    end
    old_message, v[factor] = v[factor], message
    set(v, v / old_message * message)
end

function update_value(v::Variable, factor::Node; pi=0, tau=0, value::Union{Nothing, Gaussian}=nothing)
    if value === nothing
        value = Gaussian(pi, tau)
    end
    old_message = v[factor]
    v[factor] = value * old_message / v
    return set(v, value)
end

function show(io::IO, v::Variable)
    print(io, "< Variable $(v.val) with $(length(v.messages)) connections >")
end

function show(io::IO, v::Node)
    print(io, "< $(typeof(v))  with $(length(v.vars)) connections >")
end

struct Factor <: Node
    vars::Vector{Variable}
    function Factor(variables::Vector{Variable})
        r = new(variables)
        for var in variables
            var[r] = Gaussian()
        end
        return r
    end
end

@inline down(f::Factor) = 0.0
@inline up(f::Factor) = 0.0

function var(f::Node) 
    @assert length(f.vars) == 1
    return f.vars[1]
end


struct PriorFactor <: Node
    vars::Vector{Variable}
    val::Gaussian
    dynamic::Float64
    function PriorFactor(var::Variable, val::GaussianItem, dynamic::Float64=0.0)
        r = new([var], Gaussian(val.pi, val.tau), dynamic)
        var[r] = Gaussian()
        return r
    end
end

function down(f::PriorFactor)
    sig = sqrt(sigma(f.val)^2 + f.dynamic^2)
    value = Gaussian(; mu=mu(f.val), sigma=sig)
    #println("update_value from down(::PriorFactor): $(var(f)), $f; value=$value")
    return update_value(var(f), f; value=value)
end


struct LikelihoodFactor <: Node
    vars::Vector{Variable}
    mean::Variable
    value::Variable
    variance::Float64
    function LikelihoodFactor(mean_var::Variable, value_var::Variable, variance::Float64)
        r = new([mean_var, value_var], mean_var, value_var, variance)
        mean_var[r] = Gaussian()
        value_var[r] = Gaussian()
        return r
    end
end

@inline function calc_a(f::LikelihoodFactor, g::Gaussian)
    return 1 ./ (1. + f.variance * g.pi)
end

@inline function calc_a(f::LikelihoodFactor, var::Variable)
    return 1. / (1. + f.variance * var.val.pi)
end

function down(f::LikelihoodFactor)
    # update value
    msg = f.mean / f.mean[f]
    a = calc_a(f, msg)
    #println("update_message from down(f::LikelihoodFactor): $(f.value), $f, pi = $(a * msg.pi), tau = $(a * msg.tau)")    
    return update_message(f.value, f; pi = a * msg.pi, tau = a * msg.tau)
end

function up(f::LikelihoodFactor)
    # update mean
    msg = f.value / f.value[f]
    a = calc_a(f, msg)
    #println("update_message from up(f::LikelihoodFactor): $(f.mean), $f, pi = $(a * msg.pi), tau = $(a * msg.tau)")
    return update_message(f.mean, f; pi = a * msg.pi, tau = a * msg.tau)
end


struct SumFactor <: Node
    vars::Vector{Variable}
    sum::Variable
    terms::Vector{Variable}
    coeffs::Vector{Float64}
    function SumFactor(sum_var::Variable, term_vars::Vector{Variable}, coeffs::Vector{T}) where T <: Real
        r = new(vcat([sum_var], term_vars), sum_var, term_vars, coeffs)
        sum_var[r] = Gaussian()
        for var in term_vars
            var[r] = Gaussian()
        end
        return r 
    end
end

function update(f::SumFactor, var, vals, msgs, coeffs)
    pi_inv = 0.
    _mu = 0.
    for (val, msg, coeff) in zip(vals, msgs, coeffs)
        div = val / msg
        _mu += coeff * mu(div)
        if pi_inv == Inf
            continue
        end
        if div.pi == 0
            pi_inv = Inf
        else
            pi_inv += coeff ^ 2 / div.pi
        end
    end
    
    pi = 1. / pi_inv
    tau = pi * _mu
    #println("update_message from update(::SumFactor, ...): $var, $f; pi=$pi, tau=$tau")
    return update_message(var, f; pi=pi, tau=tau)
end

function down(f::SumFactor)
    vals = f.terms
    msgs = [var[f] for var in vals]
    #println("update from down(::SumFactor): $f, $(f.sum), $vals, $msgs, $(f.coeffs)")
    return update(f, f.sum, vals, msgs, f.coeffs)
end

function up(f::SumFactor; index=1)
    coeff = f.coeffs[index]
    coeffs = Array{Float64}(undef, length(f.coeffs))
    for (i, c)  in enumerate(f.coeffs)
        if coeff == 0
            coeffs[i] = 0
        elseif i == index
            coeffs[i] = 1. / coeff
        else
            coeffs[i] = -c / coeff
        end
    end
    vals = copy(f.terms)
    vals[index] = f.sum
    msgs = [var[f] for var in vals]
    #println("update from up(::SumFactor): $f, $(f.terms[index]), $vals, $msgs, $coeffs")
    return update(f, f.terms[index], vals, msgs, coeffs)
end


struct TruncateFactor <: Node
    vars::Vector{Variable}
    v_func::Function
    w_func::Function
    draw_margin::Float64
    function TruncateFactor(var, v_func, w_func, draw_margin)
        r = new([var], v_func, w_func, draw_margin)
        var[r] = Gaussian()
        return r 
    end
end

function up(f::TruncateFactor)
    val = var(f)
    msg = var(f)[f]
    div = val / msg
    sqrt_pi = sqrt(div.pi)
    args = (div.tau / sqrt_pi, f.draw_margin * sqrt_pi)
    v = f.v_func(args...)
    w = f.w_func(args...)
    denom = (1. - w)
    pi, tau = div.pi / denom, (div.tau + sqrt_pi * v) / denom
    #println("update_value from up(::TruncateFactor): $val, $f; pi=$pi, tau=$tau")
    return update_value(val, f; pi=pi, tau=tau)
end

