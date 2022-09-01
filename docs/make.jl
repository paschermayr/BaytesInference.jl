using BaytesInference
using Documenter

DocMeta.setdocmeta!(BaytesInference, :DocTestSetup, :(using BaytesInference); recursive=true)

makedocs(;
    modules=[BaytesInference],
    authors="Patrick Aschermayr <p.aschermayr@gmail.com>",
    repo="https://github.com/paschermayr/BaytesInference.jl/blob/{commit}{path}#{line}",
    sitename="BaytesInference.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://paschermayr.github.io/BaytesInference.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Introduction" => "intro.md",
    ],
)

deploydocs(;
    repo="github.com/paschermayr/BaytesInference.jl",
    devbranch="main",
)
