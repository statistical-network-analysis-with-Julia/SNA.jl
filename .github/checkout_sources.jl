using TOML

# CI checks out this package as SNA.jl beside its local [sources] dependencies.
# The clone list comes from committed metadata instead of drifting YAML lists.
root = normpath(joinpath(@__DIR__, "..", ".."))
function checkout_sources(package_dir)
    project = TOML.parsefile(joinpath(package_dir, "Project.toml"))
    for (_, source) in get(project, "sources", Dict())
        haskey(source, "path") || continue
        target = normpath(joinpath(package_dir, source["path"]))
        target == package_dir && continue
        relpath(target, root) == basename(target) ||
            error("local dependency must be a sibling checkout: $target")
        if !isdir(target)
            repo = "https://github.com/statistical-network-analysis-with-Julia/$(basename(target)).git"
            run(`git clone --depth 1 --quiet $repo $target`)
            checkout_sources(target)
        end
    end
end
checkout_sources(joinpath(root, "SNA.jl"))
