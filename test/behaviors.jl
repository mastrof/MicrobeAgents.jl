using MicrobeAgents, Test
using MicrobeAgents: _foreach, _prod, copybehavior, copybehaviors

# types and hook methods must be defined at top level, not inside @testset
struct Dummy end
mutable struct Counter
    n::Int
end

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
end
