module EffectSizes

export
    EffectSize,
    CohenD,
    HedgeG,
    effectsize,
    confint,
    quantile

import Statistics: mean, std
import Distributions
import Distributions: quantile
import HypothesisTests
import HypothesisTests: confint
import Base.Grisu: PRECISION

include("confint.jl")
include("effectsize.jl")

const EffectSize = CohenD

end # module
