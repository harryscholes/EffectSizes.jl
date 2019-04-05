"""
    AbstractEffectSize

An abstract type to represent an effect size.

Effect | Effect size
:---|---:
Small | 0.2
Medium | 0.5
Large | 0.8

Subtypes implement:

Method | Description
:--- | :---
`effectsize` | returns the effect size index
`confint` | returns the confidence interval
`quantile` | returns the confidence interval quantile
"""
abstract type AbstractEffectSize end

HypothesisTests.confint(es::AbstractEffectSize) = es.ci
Distributions.quantile(es::AbstractEffectSize) = quantile(confint(es))

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

function correction(n::Integer)
    n > 1 || throw(DomainError(n))
    ((n-3)/(n-2.25)) * √((n-2)/n)
end

effectsize(mx::T, my::T, s::T) where T<:Real = ((mx-my) / s)
effectsize(mx::T, my::T, s::T, n::Integer) where T<:Real = ((mx-my) / s) * correction(n)

"""
    CohenD(xs, ys[, bootstrap]; [quantile=0.95])

Calculate Cohen's ``d`` effect size index between two vectors `xs` and `ys`.

A confidence interval for the effect size is calculated at the `quantile` quantile. If
`bootstrap` is provided, the confidence interval is calculated by resampling from `xs`
and `ys` `bootstrap` times.

```math
    d = \\frac{m_A - m_B}{s}
```

where ``m`` is the mean and ``s`` is the pooled standard deviation:

```math
    s = \\sqrt{\\frac{(n_A - 1) s_A^2 + (n_B - 1) s_B^2}{n_A + n_B - 2}}
```

If ``m_A`` > ``m_B``, ``d`` will be positive and if ``m_A`` < ``m_B``, ``d`` will be negative.

!!! note

    `HedgeG` outperforms `CohenD` when sample sizes are < 20.

# Examples

```julia
xs = randn(100000)
ys = randn(100000) .+ 0.01

using EffectSizes
CohenD(xs, ys)

using HypothesisTests
EqualVarianceTTest(xs, ys)
```
"""
struct CohenD{T<:Real} <: AbstractEffectSize
    d::T
    ci::ConfidenceInterval{T}
end

function CohenD(xs::AbstractVector{T}, ys::AbstractVector{T};
                quantile::Float64=0.95) where T<:Real
    d = cohend(xs, ys)
    ci = ConfidenceInterval(xs, ys, d; quantile=quantile)
    CohenD(d, ci)
end

function CohenD(xs::AbstractVector{T}, ys::AbstractVector{T},
                bootstrap::Integer; quantile::Float64=0.95) where T<:Real
    d = cohend(xs, ys)
    ci = ConfidenceInterval(xs, ys, cohend, bootstrap; quantile=quantile)
    CohenD(d, ci)
end

cohend(xs, ys) = effectsize(mean(xs), mean(ys), pooledstd1(xs, ys), length(xs)+length(ys))

effectsize(es::CohenD) = es.d

"""
    HedgeG(xs, ys[, bootstrap]; [quantile=0.95])

Calculate Hedge's ``g`` effect size index between two vectors `xs` and `ys`.

A confidence interval for the effect size is calculated at the `quantile` quantile. If
`bootstrap` is provided, the confidence interval is calculated by resampling from `xs`
and `ys` `bootstrap` times.

```math
    g = \\frac{m_A - m_B}{s}
```

where ``m`` is the mean and ``s`` is the pooled standard deviation:

```math
    s = \\sqrt{\\frac{(n_A - 1) s_A^2 + (n_B - 1) s_B^2}{n_A + n_B}}
```

If ``m_A`` > ``m_B``, ``g`` will be positive and if ``m_A`` < ``m_B``, ``g`` will be
negative.

!!! note

    `HedgeG` outperforms `CohenD` when sample sizes are < 20.
"""
struct HedgeG{T<:Real} <: AbstractEffectSize
    g::T
    ci::ConfidenceInterval{T}
end

function HedgeG(xs::AbstractVector{T}, ys::AbstractVector{T};
                quantile::Float64=0.95) where T<:Real
    g = hedgeg(xs, ys)
    ci = ConfidenceInterval(xs, ys, g; quantile=quantile)
    HedgeG(g, ci)
end

function HedgeG(xs::AbstractVector{T}, ys::AbstractVector{T},
                bootstrap::Integer; quantile::Float64=0.95) where T<:Real
    g = hedgeg(xs, ys)
    ci = ConfidenceInterval(xs, ys, hedgeg, bootstrap; quantile=quantile)
    HedgeG(g, ci)
end

hedgeg(xs, ys) = effectsize(mean(xs), mean(ys), pooledstd2(xs, ys), length(xs)+length(ys))

effectsize(es::HedgeG) = es.g

"""
    GlassΔ(treatment, control[, bootstrap]; [quantile=0.95])

Calculate Glass's ``Δ`` effect size index between two vectors `treatment` and `control`.

A confidence interval for the effect size is calculated at the `quantile` quantile. If
`bootstrap` is provided, the confidence interval is calculated by resampling from `xs`
and `ys` `bootstrap` times.

```math
    Δ = \\frac{m_T - m_C}{s_C}
```

where ``m`` is the mean, ``s`` is the standard deviation, ``T`` is the treatment group and
``C`` is the control group.

If ``m_T`` > ``m_C``, ``Δ`` will be positive and if ``m_T`` < ``m_C``, ``Δ`` will be negative.

!!! note

    `GlassΔ` should be used when the standard deviations between the two groups are very
    different.
"""
struct GlassΔ{T<:Real} <: AbstractEffectSize
    Δ::T
    ci::ConfidenceInterval{T}
end

function GlassΔ(xs::AbstractVector{T}, ys::AbstractVector{T};
                quantile::Float64=0.95) where T<:Real
    Δ = glassΔ(xs, ys)
    ci = ConfidenceInterval(xs, ys, Δ; quantile=quantile)
    GlassΔ(Δ, ci)
end

function GlassΔ(xs::AbstractVector{T}, ys::AbstractVector{T},
                bootstrap::Integer; quantile::Float64=0.95) where T<:Real
    Δ = glassΔ(xs, ys)
    ci = ConfidenceInterval(xs, ys, glassΔ, bootstrap; quantile=quantile)
    GlassΔ(Δ, ci)
end

glassΔ(xs, ys) = effectsize(mean(xs), mean(ys), std(ys))

effectsize(es::GlassΔ) = es.Δ
