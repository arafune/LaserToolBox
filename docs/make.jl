using LaserToolBox
using Documenter

DocMeta.setdocmeta!(LaserToolBox, :DocTestSetup, :(using LaserToolBox); recursive = true)

makedocs(;
    modules = [LaserToolBox],
    authors = "Ryuichi ARafune",
    sitename = "LaserToolBox.jl",
    format = Documenter.HTML(; edit_link = "main", assets = String[]),
    pages = ["Home" => "index.md"],
    warnonly = [:missing_docs, :cross_references],
)
