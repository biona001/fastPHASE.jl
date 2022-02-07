using Documenter
using fastPHASE

makedocs(
    sitename = "fastPHASE",
    format = Documenter.HTML(),
    modules = [fastPHASE],
    pages = [
        "Home" => "index.md",
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo   = "github.com/biona001/fastPHASE.jl.git",
    target = "build"
)