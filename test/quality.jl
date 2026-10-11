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

@testset "JET obstacles stepping" begin
    space = ContinuousSpace((50.0, 50.0); periodic = (true, false))
    sets = (
        Obstacles(space, [ObstacleSphere((25.0, 25.0), 5.0)]; contact = Slide()),
        Obstacles(space, [ObstacleCapsule((5.0, 5.0), (15.0, 10.0), 1.0)]; contact = Reflect()),
        Walls(space; contact = Randomize()),
    )
    model = StandardABM(Microbe{2,2,Tuple{}}, space, 0.1; properties = Dict(:obstacles => sets))
    motility = RunTumble(run_speed = [30.0], run_duration = 1.0, angle = Isotropic(2))
    m = add_agent!(SVector(10.0, 40.0), model; motility, radius = 0.5)
    # the translation path changed by obstacles (reorient_step! has known
    # runtime dispatch on motile states, out of scope here)
    JET.@test_opt target_modules = (MicrobeAgents,) move_step!(m, model)
end
