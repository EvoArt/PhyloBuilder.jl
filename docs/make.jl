using PhyloBuilder
using Documenter

DocMeta.setdocmeta!(PhyloBuilder, :DocTestSetup, :(using PhyloBuilder); recursive=true)

makedocs(;
    modules=[PhyloBuilder],
    authors="Arthur Newbury",
    sitename="PhyloBuilder.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://EvoArt.github.io/PhyloBuilder.jl",
        edit_link="master",
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => "examples.md",
        "API" => "api.md",
    ],
    # The temporary bundled PhylogenyIO module has its own separate manual.
    checkdocs=:none,
)

deploydocs(;
    repo="github.com/EvoArt/PhyloBuilder.jl",
    devbranch="master",
)
