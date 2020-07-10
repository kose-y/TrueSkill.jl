module TrueSkill
    import Base: show, isless, *, /
    export TrueSkillEnv, Rating, expose, rate, mu, sigma
    include("gaussian.jl")
    include("factorgraph.jl")
    include("distribution.jl")
    include("mainbody.jl")
end
