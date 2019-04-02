"""
    AbstractEffectSize

An abstract type to represent an effect size.

Effect | Effect size
---|---
Small | 0.2
---|---
Medium | 0.5
---|---
Large | 0.8

Subtypes implement:
    - `effectsize`, returns the effect size index
    - `confint`, returns the confidence interval
    - `quantile`, returns the confidence interval quantile
"""
abstract type AbstractEffectSize end

HypothesisTests.confint(es::AbstractEffectSize) = es.ci
Distributions.quantile(es::AbstractEffectSize) = confint(es).quantile

function Base.show(io::IO, es::AbstractEffectSize)
    print(io, round(effectsize(es), digits=PRECISION), ", ", confint(es))
end

# Used by: CohenD
function pooledstd1(xs::AbstractVector{T}, ys::AbstractVector{T}) where T<:Real
    nx = length(xs)
    ny = length(ys)
    √(((nx-1)*std(xs)^2 + (ny-1)*std(ys)^2)/(nx+ny-2))
end

# Used by: HedgeG
function pooledstd2(xs::AbstractVector{T}, ys::AbstractVector{T}) where T<:Real
    nx = length(xs)
    ny = length(ys)
    √(((nx-1)*std(xs)^2 + (ny-1)*std(ys)^2)/(nx+ny))
end

correction(n::Integer) = ((n-3)/(n-2.25)) * √((n-2)/n)

"""
    CohenD(xs, ys[, bootstrap]; [quantile=0.95])

Calculate Cohen's ``d`` effect size index between two vectors `xs` and `ys`.

```math
    d = \\frac{m_A - m_B}{s}
```

where ``m`` is the mean and ``s`` is the pooled standard deviation:

```math
    s = \\sqrt{\\frac{(n_A - 1) s_A^2 + (n_B - 1) s_B^2}{n_A + n_B - 2}}
```

If ``m_A`` > ``m_B``, ``d`` will be positive and if ``m_A`` < ``m_B``, ``d`` will be negative.

# Examples

```julia
xs = rand(100000)
ys = rand(100000) *.99

EffectSize(xs, ys)

using HypothesisTests
pvalue(EqualVarianceTTest(xs, ys))
```
"""
struct CohenD{T<:Real} <: AbstractEffectSize
    d::T
    ci::ConfidenceInterval{T}
end

function CohenD(xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real};
                quantile::Float64=0.95)
    d = cohend(xs, ys)
    ci = ConfidenceInterval(xs, ys, d; quantile=quantile)
    CohenD(d, ci)
end

function CohenD(xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real},
                bootstrap::Integer; quantile::Float64=0.95)
    d = cohend(xs, ys)
    ci = ConfidenceInterval(xs, ys, cohend, bootstrap; quantile=quantile)
    CohenD(d, ci)
end

function cohend(mx::T, my::T, s::T, n::Integer) where T<:Real
    d = ((mx-my) / s) * correction(n)
end

function cohend(xs::AbstractVector{T}, ys::AbstractVector{T}) where T<:Real
    cohend(mean(xs), mean(ys), pooledstd1(xs, ys), length(xs)+length(ys))
end

effectsize(es::CohenD) = es.d

"""
    HedgeG(xs, ys[, bootstrap]; [quantile=0.95])

Calculate Hedge's ``g`` effect size index between two vectors `xs` and `ys`.

```math
    g = \\frac{m_A - m_B}{s}
```

where ``m`` is the mean and ``s`` is the pooled standard deviation:

```math
    s = \\sqrt{\\frac{(n_A - 1) s_A^2 + (n_B - 1) s_B^2}{n_A + n_B}}
```

If ``m_A`` > ``m_B``, ``d`` will be positive and if ``m_A`` < ``m_B``, ``d`` will be negative.
"""
struct HedgeG{T<:Real} <: AbstractEffectSize
    g::T
    ci::ConfidenceInterval{T}
end

function HedgeG(xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real};
                quantile::Float64=0.95)
    g = hedgeg(xs, ys)
    ci = ConfidenceInterval(xs, ys, g; quantile=quantile)
    HedgeG(g, ci)
end

function HedgeG(xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real},
                bootstrap::Integer; quantile::Float64=0.95)
    g = hedgeg(xs, ys)
    ci = ConfidenceInterval(xs, ys, hedgeg, bootstrap; quantile=quantile)
    HedgeG(g, ci)
end

function hedgeg(xs::AbstractVector{T}, ys::AbstractVector{T}) where T<:Real
    cohend(mean(xs), mean(ys), pooledstd2(xs, ys), length(xs)+length(ys))
end

effectsize(es::HedgeG) = es.g
