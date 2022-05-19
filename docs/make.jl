using Plastic
using Documenter

DocMeta.setdocmeta!(Plastic, :DocTestSetup, :(using Plastic); recursive=true)

makedocs(;
    modules=[Plastic],
    authors="iamzhtr@hotmail.com",
    repo="https://github.com/HetaoZ/Plastic.jl/blob/{commit}{path}#{line}",
    sitename="Plastic.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://HetaoZ.github.io/Plastic.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/HetaoZ/Plastic.jl",
    devbranch="main",
)
