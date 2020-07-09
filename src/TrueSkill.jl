module TrueSkill
    import Base: show, isless, *, /
    export TrueSkillEnv, Rating, expose, rate, rate_1vs1 , rate_freeforall
    include("gaussian.jl")
    include("factorgraph.jl")
    include("distribution.jl")
    include("mainbody.jl")
end
