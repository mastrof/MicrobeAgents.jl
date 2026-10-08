using Aqua, JET

@testset "Aqua" begin
    Aqua.test_all(MicrobeAgents;
        # `add_agent!(pos, Type, model, ...)` is inherently ambiguous with Agents' catch-all methods
        ambiguities = (exclude = [MicrobeAgents.Agents.add_agent!],),
        # `position` is extended on `SVector`/`NTuple` on purpose
        piracies = (treat_as_own = [Base.position],),
    )
end

@testset "JET" begin
    rep = JET.report_package(MicrobeAgents;
        target_defined_modules = true, toplevel_logger = nothing
    )
    # known false positives from macro-generated code:
    # - `@agent`+`@kwdef` constructors without `D` (defaults referring to `D`)
    # - LightSumTypes' generic `copy` over the motile states
    known(r) = let msg = sprint(JET.print_report, r)
        occursin(r"MicrobeAgents\.D` is not defined", msg) ||
        occursin(r"no matching method found `copy\(::(\S*\.)?(Run|Turn)State\)`", msg)
    end
    reports = filter(!known, JET.get_reports(rep))
    isempty(reports) || foreach(r -> println(sprint(JET.print_report, r)), reports)
    @test isempty(reports)
end
