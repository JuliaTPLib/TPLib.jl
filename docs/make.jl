using TPLib, Documenter

makedocs(
    modules = [TPLib],
    format = Documenter.HTML(; assets = ["assets/custom.css"]),
    sitename = "TPLib.jl",
    pages = [
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Usage" => "usage.md",
        "Reference" => "reference.md"
    ]
)

deploydocs(
    repo = "github.com/JuliaTPLib/TPLib.jl.git"
)
