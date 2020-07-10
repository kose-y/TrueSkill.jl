import SpecialFunctions: erfc, erfcinv

"""
inverse of cdf.
"""
@inline function ppf(x)
    return -  sqrt(2) * erfcinv(2 * x)
end

@inline function pdf(x)
    return (1. / (sqrt(2 * pi) ) *
        exp(-(((x)) ^ 2 / 2.0)))
end

@inline function cdf(x)
    return 0.5 * erfc(-(x) / (sqrt(2)))
end
