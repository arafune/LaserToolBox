using LaserToolBox
using Coverage
using Test
if Base.JLOptions().code_coverage == 0
    println("Restarting with --code-coverage for coverage measurement...")
    run(
        `$(Base.julia_cmd()) --project=$(Base.active_project()) --code-coverage test/runtests_local.jl`,
    )
    exit()
end

include("runtests.jl")

results = process_folder("src")
LCOV.writefile("lcov.info", results)
run(`genhtml  --ignore-errors range  -o coverage_html lcov.info`)
