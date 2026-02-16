using Documenter
using SNA

DocMeta.setdocmeta!(SNA, :DocTestSetup, :(using SNA); recursive=true)

makedocs(
    sitename = "SNA.jl",
    modules = [SNA],
    authors = "Simone Santoni",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://Statistical-network-analysis-with-Julia.github.io/SNA.jl",
        edit_link = "main",
    ),
    repo = "https://github.com/Statistical-network-analysis-with-Julia/SNA.jl/blob/{commit}{path}#{line}",
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "User Guide" => [
            "Centrality" => "guide/centrality.md",
            "Network Measures" => "guide/measures.md",
            "Cohesion" => "guide/cohesion.md",
            "Structural Equivalence" => "guide/equivalence.md",
        ],
        "API Reference" => [
            "Centrality" => "api/centrality.md",
            "Measures" => "api/measures.md",
            "Cohesion" => "api/cohesion.md",
            "Equivalence" => "api/equivalence.md",
        ],
    ],
    warnonly = [:missing_docs, :docs_block],
)

deploydocs(
    repo = "github.com/Statistical-network-analysis-with-Julia/SNA.jl.git",
    devbranch = "main",
    versions = [
        "stable" => "dev", # serve dev docs at /stable until a release is tagged
        "dev" => "dev",
    ],
    push_preview = true,
)
