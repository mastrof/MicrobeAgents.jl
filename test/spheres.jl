using MicrobeAgents, Test
using Random

@testset "Spheres" begin
    L = 100.0
    dt = 1.0
    for D in 1:3
        extent = ntuple(_ -> L, D)
        model = ABM(Microbe{D}, extent, dt)

        c1 = extent ./ 2
        r1 = 30.0
        s1 = HyperSphere(c1, r1)
        @test Tuple(s1.center) == c1
        @test s1.r == r1

        r2 = 5.0
        c2 = ntuple(i -> i==1 ? L-10 : L/2, D)
        s2 = HyperSphere(c2, r2)
        r3 = 20.0
        c3 = ntuple(i -> i==1 ? L-10 : L/2, D)
        s3 = HyperSphere(c3, r3)
        @test distance(s1,s2,model) ≈ L/2-10
        @test distance(s1,s3,model) ≈ L/2-10
        @test distance(s2,s3,model) ≈ 0
        @test ~contact(s1,s2,model)
        @test contact(s1,s3,model) 
        @test contact(s2,s3,model)

        x1 = c1 .+ ntuple(i -> i==D ? r1+1 : 0, D)
        x2 = c1 .+ ntuple(i -> i==D ? r1-1 : 0, D)
        v = ntuple(i -> i==D ? 1.0 : 0.0, D)
        add_agent!(x1, model; vel=.-v)
        add_agent!(x2, model; vel=.+v)

        @test ~contact(model[1],s1,model)
        @test contact(model[2],s1,model)
        @test_throws KeyError is_encounter(model[1],s1,model)
        # add test with proper encounter checking
    end
end