using EffectSizes
using Test

@testset "EffectSizes.jl" begin
    include("test_confint.jl")
    include("test_effectsize.jl")
end
