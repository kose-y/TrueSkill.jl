abstract type GaussianItem end

mutable struct Gaussian <: GaussianItem
    pi::Float64
    tau::Float64
end

function Gaussian(; mu::Float64=0.0, sigma::Float64=0.0)
    if sigma == 0
        return Gaussian(0.0, 0.0)
    end
    pi = sigma^(-2)
    tau = pi * mu
    Gaussian(pi, tau)
end

@inline function mu(g::Gaussian)
    if g.pi != 0
        return g.tau / g.pi
    else
        return 0
    end
end

@inline function sigma(g::Gaussian)
    if g.pi != 0
        return(sqrt(1/g.pi))
    else
        return Inf
    end
end

@inline function isless(a::Gaussian, b::Gaussian)
    isless(mu(a), mu(b))
end

@inline (*)(a::Gaussian, b::Gaussian) = Gaussian(a.pi + b.pi, a.tau + b.tau)
@inline (/)(a::Gaussian, b::Gaussian) = Gaussian(a.pi - b.pi, a.tau - b.tau)

function show(io::IO, g::Gaussian)
    print(io, "Gaussian(mu=$(mu(g)), sigma=$(sigma(g)))")
end