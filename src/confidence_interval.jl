bootstrapsample(xs::AbstractVector) = @views xs[rand(1:length(xs), length(xs))]

function twotailedquantile(quantile::AbstractFloat)
    0. ≤ quantile ≤ 1. || throw(DomainError(quantile))
    lq = (1 - quantile) / 2
    uq = quantile + lq
    return lq, uq
end

"""
    AbstractConfidenceInterval{T<:Real}

A type representing a confidence interval.

Subtypes implement:

Method | Description
:--- | :---
`lower` | returns the lower bound
`upper` | returns the upper bound
`quantile` | returns the quantile
"""
abstract type AbstractConfidenceInterval{T<:Real} end

"""
    lower(ci::AbstractConfidenceInterval{T}) -> T

Return the lower bound of a confidence interval.
"""
lower

"""
    upper(ci::AbstractConfidenceInterval{T}) -> T

Return the upper bound of a confidence interval.
"""
upper

"""
    quantile(ci::AbstractConfidenceInterval) -> Float64

Returns the quantile of a confidence interval.
"""
quantile

"""
    confint(ci::AbstractConfidenceInterval{T}) -> Tuple{T,T}

Return the lower and upper bounds of a confidence interval.
"""
confint(ci::AbstractConfidenceInterval) = (lower(ci), upper(ci))

"""
    ConfidenceInterval(lower, upper, quantile)

A type representing the `lower` and `upper` bounds of an effect size confidence
interval at a specified `quantile`.

    ConfidenceInterval(xs, ys, es; quantile)

Calculate a confidence interval for the effect size `es` between two vectors `xs` and `ys`
at a specified `quantile`.
"""
struct ConfidenceInterval{T<:Real} <: AbstractConfidenceInterval{T}
    lower::T
    upper::T
    quantile::Float64

    function ConfidenceInterval(l::T, u::T, q::Float64) where T<:Real
        l ≤ u || throw(ArgumentError("l > u"))
        0. ≤ q ≤ 1. || throw(DomainError(q))
        new{T}(l, u, q)
    end
end

function ConfidenceInterval(
    xs::AbstractVector{T},
    ys::AbstractVector{T},
    es::T;
    quantile::AbstractFloat,
) where T<:Real

    0. ≤ quantile ≤ 1. || throw(DomainError(quantile))
    nx = length(xs)
    ny = length(ys)
    σ² = √(((nx + ny) / (nx * ny)) + (es^2 / 2(nx + ny)))
    _, uq = twotailedquantile(quantile)
    z = Distributions.quantile(Normal(), uq)
    ci = z * √σ²

    return ConfidenceInterval(
        es - ci,
        es + ci,
        quantile,
    )
end

lower(ci::ConfidenceInterval) = ci.lower
upper(ci::ConfidenceInterval) = ci.upper
quantile(ci::ConfidenceInterval) = ci.quantile

"""
    BootstrapConfidenceInterval(lower, upper, quantile, bootstrap)

A type representing the `lower` and `upper` bounds of an effect size confidence
interval at a specified `quantile` with `bootstrap` resamples.

    BootstrapConfidenceInterval(f, xs, ys, bootstrap; quantile)

Calculate a bootstrap confidence interval between two vectors `xs` and `ys` at a specified
`quantile` by applying `f` to `bootstrap` resamples of `xs` and `ys`.
"""
struct BootstrapConfidenceInterval{T<:Real} <: AbstractConfidenceInterval{T}
    lower::T
    upper::T
    quantile::Float64
    bootstrap::Int64

    function BootstrapConfidenceInterval(l::T, u::T, q::Float64, b::Int64) where T<:Real
        l ≤ u || throw(ArgumentError("l > u"))
        0. ≤ q ≤ 1. || throw(DomainError(q))
        b > 1 || throw(DomainError(b))
        new{T}(l, u, q, b)
    end
end

function BootstrapConfidenceInterval(
    f::Function,
    xs::AbstractVector{T},
    ys::AbstractVector{T},
    bootstrap::Integer = 1000;
    quantile::AbstractFloat,
    ) where T<:Real

    0. ≤ quantile ≤ 1. || throw(DomainError(quantile))
    bootstrap > 1 || throw(DomainError(bootstrap))
    es = map(_->f(bootstrapsample(xs), bootstrapsample(ys)), 1:bootstrap)
    lq, uq = twotailedquantile(quantile)

    return BootstrapConfidenceInterval(
        Distributions.quantile(es, lq),
        Distributions.quantile(es, uq),
        quantile,
        bootstrap,
    )
end

lower(ci::BootstrapConfidenceInterval) = ci.lower
upper(ci::BootstrapConfidenceInterval) = ci.upper
quantile(ci::BootstrapConfidenceInterval) = ci.quantile
