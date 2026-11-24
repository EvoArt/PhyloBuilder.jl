@testset "distances" begin
    aln = readalignment(fixture("small.fa"))

    @testset "parity with ape::dist.dna ($m)" for (m, γ, file) in (
            (:raw, nothing, "raw"), (:p, nothing, "p"), (:JC69, nothing, "JC69"),
            (:K80, nothing, "K80"), (:F81, nothing, "F81"), (:TN93, nothing, "TN93"),
            (:logdet, nothing, "logdet"), (:paralin, nothing, "paralin"),
            (:JC69, 0.5, "JC69_g05"), (:K80, 0.5, "K80_g05"), (:F81, 0.7, "F81_g07"),
            (:TN93, 0.7, "TN93_g07"))
        d = distances(aln; model=m, gamma=γ, progress=false)
        @test d.d ≈ readmatrix("dist_$file.txt") atol = 1e-12
        @test d.model === m
    end

    @testset "DistanceMatrix" begin
        d = distances(aln; progress=false)
        @test d isa AbstractMatrix{Float64}
        @test size(d) == (8, 8) && ntaxa(d) == 8 && taxa(d) == taxa(aln)
        @test d == d' && all(iszero, d[i, i] for i in 1:8)
        @test d["s1", "s3"] == d[1, 3]
        @test_throws KeyError d["s1", "nope"]
        @test copy(d) == d && copy(d).names == d.names
        @test_throws ArgumentError DistanceMatrix(["a", "b"], zeros(3, 3))
        @test_throws ArgumentError DistanceMatrix(["a", "a"], zeros(2, 2))
    end

    @testset "threads give the same result" begin
        @test distances(aln; model=:TN93, threaded=true, progress=false) ==
              distances(aln; model=:TN93, threaded=false, progress=false)
        rng = MersenneTwister(6)                     # more than 20 taxa: the threaded loop
        many = Alignment(["t$i" for i in 1:30], [join(rand(rng, "ACGT", 60)) for _ in 1:30])
        @test distances(many; threaded=true, progress=false) == distances(many; threaded=false, progress=false)
    end

    @testset "edge cases" begin
        @test_throws ArgumentError distances(aln; model=:HKY, progress=false)
        @test_throws ArgumentError distances(aln; model=:logdet, gamma=1.0, progress=false)
        # no shared sites → NaN; saturated → Inf
        a = Alignment(["x", "y", "z"], ["AC--", "--GT", "TGCA"])
        d = distances(a; progress=false)
        @test isnan(d["x", "y"])
        @test d["x", "z"] == Inf                    # p = 1 > 3/4
        @test distances(a; model=:p, progress=false)["x", "z"] == 1.0
        # identical sequences
        b = Alignment(["x", "y"], ["ACGTACGT", "ACGTACGT"])
        for m in (:JC69, :K80, :F81, :TN93, :logdet, :paralin)
            @test distances(b; model=m, progress=false)[1, 2] ≈ 0 atol = 1e-12
        end
    end

    @testset "patristic" begin
        t = parsenewick("((a:1,b:2):0.5,(c:1,d:1):1);")
        p = patristic(t)
        @test p["a", "b"] == 3 && p["a", "c"] == 3.5 && p["c", "d"] == 2
        @test p.model === :patristic
        # without branch lengths every branch counts 1
        @test patristic(parsenewick("((a,b),c);"))["a", "c"] == 3
    end
end
