bootstrapsample(xs::AbstractVector) = xs[rand(1:length(xs), length(xs))]

function twotailedquantile(quantile::AbstractFloat)
    lq = (1-quantile)/2
    uq = quantile+lq
    lq, uq
end

"""
    ConfidenceInterval(lower, upper, quantile)

A type representing the `lower` lower and `upper` upper bounds of an effect size confidence
interval at the `qualtile` quantile.

    ConfidenceInterval(xs, ys, d; quantile)
    ConfidenceInterval(xs, ys; quantile, bootstrap)

Calculates the effect size confidence interval between two vectors `xs` and `ys` at the
`quantile` quantile.

If `d` is provided, the confidence interval is calculated using the standard normal
distribution. If `bootstrap` is provided, the confidence interval is calculated by
resampling from `xs` and `ys` `bootstrap` times.
"""
struct ConfidenceInterval{T<:Real}
    lower::T
    upper::T
    quantile::Float64
end

function ConfidenceInterval(xs::AbstractVector{T}, ys::AbstractVector{T}, es::T;
                            quantile::AbstractFloat) where T<:Real
    0. ≤ quantile ≤ 1. || throw(DomainError(quantile))
    nx = length(xs)
    ny = length(ys)
    σ² = √(((nx+ny)/(nx*ny)) + (es^2 / 2(nx+ny)))
    _, uq = twotailedquantile(quantile)
    z = Distributions.quantile(Normal(), uq)
    ci = z*√σ²
    ConfidenceInterval(es-ci, es+ci, quantile)
end

function ConfidenceInterval(xs::AbstractVector{T}, ys::AbstractVector{T}, f::Function,
                            bootstrap::Integer=1000; quantile::AbstractFloat) where T<:Real
    0. ≤ quantile ≤ 1. || throw(DomainError(quantile))
    bootstrap > 1 || throw(DomainError(bootstrap))
    es = map(_->f(bootstrapsample(xs), bootstrapsample(ys)), 1:bootstrap)
    lq, uq = twotailedquantile(quantile)
    ConfidenceInterval(Distributions.quantile(es, lq), Distributions.quantile(es, uq),
                       quantile)
end

HypothesisTests.confint(ci::ConfidenceInterval) = ci.lower, ci.upper
Distributions.quantile(ci::ConfidenceInterval) = ci.quantile

function Base.show(io::IO, ci::ConfidenceInterval)
    print(io, ci.quantile, "CI (", round(ci.lower, digits=PRECISION), ", ",
          round(ci.upper, digits=PRECISION), ")")
end
