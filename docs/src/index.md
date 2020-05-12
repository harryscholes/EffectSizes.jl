# EffectSizes.jl

EffectSizes.jl is a Julia package for effect size measures. Confidence intervals are
assigned to effect sizes using either the normal distribution or by bootstrap resampling.
The package implements types for the following measures:

**Measure** | **Type**
---|---
Cohen's *d* | `CohenD`
Hedge's *g* | `HedgeG`
Glass's *Δ* | `GlassΔ`

## Installation

```julia
] add https://github.com/harryscholes/EffectSizes.jl
```

## Examples

```julia
julia> using Random, EffectSizes; Random.seed!(1);

julia> xs = randn(10^3);

julia> ys = randn(10^3) .+ 0.5;

julia> es = CohenD(xs, ys, quantile=0.95); # normal CI (idealised distribution)

julia> typeof(es)
CohenD{Float64,ConfidenceInterval{Float64}}

julia> effectsize(es)
-0.507…

julia> quantile(es)
0.95

julia> ci = confint(es);

julia> typeof(ci)
ConfidenceInterval{Float64}

julia> confint(ci)
(-0.924…, -0.0889…)

julia> es = CohenD(xs, ys, 10^4, quantile=0.95); # bootstrap CI (empirical distribution)

julia> effectsize(es) # effect size is the same
-0.507…

julia> typeof(es)
CohenD{Float64,BootstrapConfidenceInterval{Float64}}

julia> ci = confint(es); # confidence interval is different

julia> lower(ci)
-0.597…

julia> upper(ci)
-0.418…
```

## Index

```@index
```

## API

```@docs
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
```
