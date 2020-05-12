"""
    EffectSizes

A Julia package for effect size measures.
"""
module EffectSizes

export
    AbstractEffectSize,
    EffectSize,
    CohenD,
    HedgeG,
    GlassΔ,
    effectsize,
    confint,
    quantile,
    AbstractConfidenceInterval,
    ConfidenceInterval,
    BootstrapConfidenceInterval,
    lower,
    upper


using Statistics
using Distributions
using HypothesisTests

import Distributions: quantile
import HypothesisTests: confint

include("confidence_interval.jl")
include("effect_size.jl")

end # module
