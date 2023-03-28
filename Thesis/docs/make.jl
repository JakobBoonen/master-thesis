using Documenter, Thesis

makedocs(
    modules = [Thesis],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Jakob",
    sitename = "Thesis.jl",
    pages = Any["index.md"]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

deploydocs(
    repo = "github.com/Jakob/Thesis.jl.git",
    push_preview = true
)
