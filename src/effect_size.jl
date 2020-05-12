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
"""
abstract type AbstractEffectSize end

"""
    effectsize(es::AbstractEffectSize)

Return the effect size index.
"""
effectsize(::T) where T<:AbstractEffectSize =
    throw(ArgumentError("`effectsize` is not implemented for $(string(T))"))

"""
    confint(es::AbstractEffectSize) -> ConfidenceInterval

Return the confidence interval of an effect size as a `ConfidenceInterval` object.
"""
confint(::T) where T<:AbstractEffectSize =
    throw(ArgumentError("`confint` is not implemented for $(string(T))")) 

"""
    quantile(es::AbstractEffectSize) -> Float64

Returns the quantile of a confidence interval.
"""
quantile(es::AbstractEffectSize) = quantile(confint(es))

# Used by: CohenD
function pooledstd1(xs::AbstractVector{T}, ys::AbstractVector{T}) where T<:Real
    nx = length(xs)
    ny = length(ys)
    return √(((nx - 1) * std(xs)^2 + (ny - 1) * std(ys)^2) / (nx + ny - 2))
end

# Used by: HedgeG
function pooledstd2(xs::AbstractVector{T}, ys::AbstractVector{T}) where T<:Real
    nx = length(xs)
    ny = length(ys)
    return √(((nx - 1) * std(xs)^2 + (ny - 1) * std(ys)^2) / (nx + ny))
end

function correction(n::Integer)
    n > 1 || throw(DomainError(n))
    return ((n - 3) / (n - 2.25)) * √((n - 2) / n)
end

_effectsize(μx::T, μy::T, σ::T) where T<:Real = ((μx - μy) / σ)

_effectsize(μx::T, μy::T, σ::T, n::Integer) where T<:Real = ((μx - μy) / σ) * correction(n)

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
struct CohenD{T<:Real,CI<:AbstractConfidenceInterval{T}} <: AbstractEffectSize
    d::T
    ci::CI
end

function cohend(xs::AbstractVector{T}, ys::AbstractVector{T}) where T<:Real
    return _effectsize(
        mean(xs),
        mean(ys),
        pooledstd1(xs, ys),
        length(xs) + length(ys),
    )
end

effectsize(es::CohenD) = es.d
confint(es::CohenD) = es.ci

"""
    const EffectSize = CohenD

See [`CohenD`](@ref).
"""
const EffectSize = CohenD

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
struct HedgeG{T<:Real,CI<:AbstractConfidenceInterval{T}} <: AbstractEffectSize
    g::T
    ci::CI
end

function hedgeg(xs::AbstractVector{T}, ys::AbstractVector{T}) where T<:Real
    return _effectsize(
        mean(xs),
        mean(ys),
        pooledstd2(xs, ys),
        length(xs) + length(ys),
    )
end

effectsize(es::HedgeG) = es.g
confint(es::HedgeG) = es.ci

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
struct GlassΔ{T<:Real,CI<:AbstractConfidenceInterval{T}} <: AbstractEffectSize
    Δ::T
    ci::CI
end

function glassΔ(xs::AbstractVector{T}, ys::AbstractVector{T}) where T<:Real
    return _effectsize(
        mean(xs),
        mean(ys),
        std(ys),
    )
end

effectsize(es::GlassΔ) = es.Δ
confint(es::GlassΔ) = es.ci


# constructors

for (T, f) = [(:CohenD, cohend), (:GlassΔ, glassΔ), (:HedgeG, hedgeg)]
    @eval begin
        # Normal CI
        function $T(
            xs::AbstractVector{T},
            ys::AbstractVector{T};
            quantile::Float64=0.95,
        ) where T<:Real
            es = $f(xs, ys)
            ci = ConfidenceInterval(xs, ys, es; quantile=quantile)
            $T(es, ci)
        end

        # Bootstrap CI
        function $T(
            xs::AbstractVector{T},
            ys::AbstractVector{T},
            bootstrap::Integer;
            quantile::Float64=0.95,
        ) where T<:Real
            es = $f(xs, ys)
            ci = BootstrapConfidenceInterval($f, xs, ys, bootstrap; quantile=quantile)
            $T(es, ci)
        end
    end
end
