using MicrobeAgents, Test
using MicrobeAgents: _foreach, _prod, copybehavior, copybehaviors

# types and hook methods must be defined at top level, not inside @testset
struct Dummy end
mutable struct Counter
    n::Int
end
using Random
mutable struct Ticker
    n::Int
end
MicrobeAgents.affect!(t::Ticker, microbe, model) = (t.n += 1; nothing)
MicrobeAgents.initialize!(t::Ticker, microbe, model) = (t.n = 100; nothing)

@testset "Behaviors" begin
    @testset "Default hooks are neutral" begin
        @test isnothing(initialize!(Dummy(), nothing, nothing))
        @test isnothing(affect!(Dummy(), nothing, nothing))
        @test bias(Dummy(), nothing, nothing) === 1.0
        @test speed_factor(Dummy(), nothing) === 1.0
        w = [0.5, 0.5]
        @test isnothing(transition_weights!(w, Dummy(), nothing, nothing))
        @test w == [0.5, 0.5]
    end

    @testset "Plain functions are affect!" begin
        hits = Ref(0)
        f(microbe, model) = (hits[] += 1; :ignored)
        @test isnothing(affect!(f, nothing, nothing))
        @test hits[] == 1
        @test bias(f, nothing, nothing) === 1.0
    end

    @testset "Behavior wrapper" begin
        b = Behavior(
            bias = (m, model) -> 2.0,
            speed_factor = m -> 3.0,
            transition_weights! = (w, m, model) -> (w[1] = 0.0),
        )
        @test bias(b, nothing, nothing) === 2.0
        @test speed_factor(b, nothing) === 3.0
        w = [0.5, 0.5]
        transition_weights!(w, b, nothing, nothing)
        @test w == [0.0, 0.5]
        @test isnothing(affect!(b, nothing, nothing))
        hits = Ref(0)
        b2 = Behavior(affect! = (m, model) -> (hits[] += 1))
        affect!(b2, nothing, nothing)
        @test hits[] == 1
        @test bias(b2, nothing, nothing) === 1.0
        @test_throws MethodError Behavior(unknown = identity)
    end

    @testset "Combinators" begin
        order = Int[]
        _foreach(x -> push!(order, x), (1, 2, 3))
        @test order == [1, 2, 3]
        @test isnothing(_foreach(identity, ()))
        @test _prod(identity, ()) === 1.0
        @test _prod(identity, (2.0, 3.0)) === 6.0
        @test (@inferred _prod(x -> x isa Int ? 2.0 : 0.5, (1, "a", 3))) === 2.0
    end

    @testset "Copying" begin
        f(m, model) = nothing
        c = Counter(0)
        bs = (c, f)
        cs = copybehaviors(bs)
        @test cs[1] !== c && cs[1].n == 0
        @test cs[2] === f
        nt = (a = Counter(1), b = Behavior(bias = (m, model) -> 2.0))
        cnt = copybehaviors(nt)
        @test cnt isa NamedTuple{(:a, :b)}
        @test cnt.a !== nt.a && cnt.b === nt.b
        @test_throws ArgumentError copybehaviors(Counter(0))
    end

    @testset "Microbe with behaviors" begin
        for D in 1:3, container in (Vector, Dict)
            space = ContinuousSpace(fill(100.0, SVector{D}))
            model = StandardABM(Microbe{D}, space, 0.1; container, rng = Xoshiro(1))
            motility = RunTumble([30.0], 1.0, Isotropic(D))
            shared = (t = Ticker(0),)
            add_agent!(model; motility, behaviors = shared)
            add_agent!(model; motility, behaviors = shared)
            add_agent!(model; motility) # no behaviors: mixed population
            @test model[1] isa Microbe{D,2}
            @test behaviors(model[3]) === ()
            @test model[1].behaviors.t !== model[2].behaviors.t # independent copies
            @test shared.t.n == 0 # original untouched
            @test model[1].behaviors.t.n == 100 # initialize! ran
            run!(model, 3)
            @test model[1].behaviors.t.n == model[2].behaviors.t.n == 103
            @test findbehavior(model[1], Ticker) === model[1].behaviors.t
            @test findbehavior(model[3], Ticker) === nothing
        end
    end

    @testset "Single behavior without tuple" begin
        space = ContinuousSpace((10.0, 10.0))
        model = StandardABM(Microbe{2}, space, 0.1)
        motility = RunTumble([30.0], 1.0, Isotropic(2))
        @test_throws ArgumentError add_agent!(model; motility,
            behaviors = Behavior(bias = (m, model) -> 2.0))
    end

    @testset "Hooks in the pipeline" begin
        space = ContinuousSpace((1000.0, 1000.0))
        dt = 0.1
        # bias: product over behaviors, used in switching_probability
        model = StandardABM(Microbe{2}, space, dt)
        motility = RunTumble([30.0], 1.0, Isotropic(2))
        b2 = Behavior(bias = (m, model) -> 2.0)
        b3 = Behavior(bias = (m, model) -> 3.0)
        add_agent!(model; motility, behaviors = (b2, b3))
        @test bias(model[1], model) == 6.0
        @test switching_probability(model[1], model) ≈ 6.0 * dt / 1.0
        # speed_factor: product, applied on top of the sampled speed
        add_agent!(model; motility, speed = 10.0,
            behaviors = (Behavior(speed_factor = m -> 2.0), Behavior(speed_factor = m -> 1.5)))
        @test model[2].speed == 10.0
        @test speed(model[2]) == 30.0
        @test velocity(model[2]) ≈ direction(model[2]) .* 30.0
        # plain function as affect!
        hits = Ref(0)
        counter(m, model) = (hits[] += 1)
        add_agent!(model; motility, behaviors = (counter,))
        run!(model, 2)
        @test hits[] == 2
    end

    @testset "transition_weights! uses a scratch copy" begin
        space = ContinuousSpace((1000.0, 1000.0))
        model = StandardABM(Microbe{2}, space, 1.0; rng = Xoshiro(3))
        # RunReverseFlick: state 3 (backward run) always goes to 4 (flick);
        # the hook redirects it to state 2 (reverse) instead
        motility = RunReverseFlick([30.0], 0.0, [30.0], 0.0)
        redirect = Behavior(transition_weights! = (w, m, model) ->
            (state(motilepattern(m)) == 3 && (w[2] = 1.0; w[4] = 0.0)))
        add_agent!(model; motility, behaviors = (redirect,))
        m = model[1]
        m.motility.current_state = 3
        MicrobeAgents.update_motilestate!(m, model)
        @test state(motilepattern(m)) == 2
        @test transition_weights(motilepattern(m), 3) == [0.0, 0.0, 0.0, 1.0] # unchanged
    end
end
