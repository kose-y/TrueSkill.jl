import SpecialFunctions: erfc, erfcinv

"""
inverse of cdf.
"""
@inline function ppf(x, mu=0, sigma=1)
    return mu - sigma * sqrt(2) * erfcinv(2 * x)
end

@inline function pdf(x, mu=0, sigma=1)
    return (1. / (sqrt(2 * pi) * abs(sigma)) *
        exp(-(((x - mu) / abs(sigma)) ^ 2 / 2.0)))
end

@inline function cdf(x, mu=0, sigma=1)
    return 0.5 * erfc(-(x - mu) / (sqrt(2) * sigma))
end
