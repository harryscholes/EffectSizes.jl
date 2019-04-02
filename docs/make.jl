using Documenter, EffectSizes

makedocs(;
    modules=[EffectSizes],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/harryscholes/EffectSizes.jl/blob/{commit}{path}#L{line}",
    sitename="EffectSizes.jl",
    authors="harryscholes",
    assets=[],
)

deploydocs(;
    repo="github.com/harryscholes/EffectSizes.jl",
)
