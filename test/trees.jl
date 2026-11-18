@testset "trees" begin
    aln = readalignment(fixture("small.fa"))
    names = taxa(aln)
    d = distances(aln; model=:K80, progress=false)

    @testset "parity with ape/phangorn ($f)" for (f, file, tol) in (
            (nj, "nj_K80.nwk", 1e-12), (upgma, "upgma_K80.nwk", 1e-12),
            (wpgma, "wpgma_K80.nwk", 1e-12),
            (bionj, "bionj_K80.nwk", 1e-6))          # ape's BIONJ is single precision
        t = f in (nj, bionj) ? f(d; nonnegative=false) : f(d)
        ref = readnewick(fixture(file))
        @test ntips(t) == 8
        @test pathdist(t, names) ≈ pathdist(ref, names) atol = tol
    end

    @testset "shapes" begin
        t = nj(d)
        @test length(childnodes(t, rootnode(t))) == 3               # unrooted
        @test nnodes(t) == 2 * 8 - 2
        @test all(i -> isroot(t, i) || branchlength(t, i) >= 0, 1:nnodes(t))
        u = upgma(d)
        @test length(childnodes(u, rootnode(u))) == 2 && nnodes(u) == 2 * 8 - 1
        h = rootdistance(u)
        @test all(≈(h[tips(u)[1]]), h[tips(u)])                      # ultrametric
        # from an alignment = from its distances
        @test writenewick(nj(aln; model=:K80, progress=false)) == writenewick(nj(d))
        @test writenewick(bionj(aln; model=:F81, gamma=0.7, progress=false)) ==
              writenewick(bionj(distances(aln; model=:F81, gamma=0.7, progress=false)))
    end

    @testset "tiny inputs" begin
        one = DistanceMatrix(["a"], zeros(1, 1))
        @test ntips(nj(one)) == 1 && ntips(upgma(one)) == 1
        two = DistanceMatrix(["a", "b"], [0 2; 2 0])
        for f in (nj, bionj, upgma)
            @test patristic(f(two))["a", "b"] ≈ 2
        end
        three = DistanceMatrix(["a", "b", "c"], [0 3 4; 3 0 5; 4 5 0])
        @test patristic(nj(three)).d ≈ three.d
        @test_throws ArgumentError nj(DistanceMatrix(String[], zeros(0, 0)))
    end

    @testset "saturated pairs" begin
        s = DistanceMatrix(["a", "b", "c", "d"], [0 1 Inf 2; 1 0 2 2; Inf 2 0 1; 2 2 1 0])
        t = @test_logs (:warn, r"saturated") nj(s)
        @test ntips(t) == 4
    end

    @testset "rooting" begin
        t = nj(d)
        P = pathdist(t, names)
        m = midpoint(t)
        @test length(childnodes(m, rootnode(m))) == 2
        @test pathdist(m, names) ≈ P
        h = rootdistance(m)
        far = maximum(P)
        @test maximum(h[tips(m)]) ≈ far / 2                         # the two farthest tips are equidistant
        @test pathdist(root(t, :midpoint), names) ≈ P
        @test pathdist(nj(d; root=:midpoint), names) ≈ P

        r = root(t, "s5")
        @test "s5" in [nodelabel(r, i) for i in childnodes(r, rootnode(r))]
        @test pathdist(r, names) ≈ P
        r2 = root(t, ["s4", "s5"])
        kids = childnodes(r2, rootnode(r2))
        @test any(k -> sort([nodelabel(r2, i) for i in tips(r2) if k in PhylogenyIO.ancestors(r2, i) || i == k]) == ["s4", "s5"], kids)
        @test pathdist(r2, names) ≈ P
        @test_throws ArgumentError root(t, :nonsense)

        u = unroot(m)
        @test length(childnodes(u, rootnode(u))) == 3
        @test pathdist(u, names) ≈ P
        @test writenewick(unroot(u)) == writenewick(u)                # already unrooted

        # support values stay with their branches (splits named by the side without "a");
        # rooting on c turns the 70 branch upside down
        s = parsenewick("((a:1,b:1)90:1,(c:1,d:1)70:2,e:1);")
        rs = root(s, "c")
        all_ = Set(["a", "b", "c", "d", "e"])
        below(t, i) = Set(nodelabel(t, j) for j in tips(t) if j == i || i in PhylogenyIO.ancestors(t, j))
        side(t, i) = (b = below(t, i); "a" in b ? setdiff(all_, b) : b)
        sup = Dict(side(rs, i) => rs.data[i][:support] for i in 1:nnodes(rs) if haskey(rs.data[i], :support))
        @test sup == Dict(Set(["c", "d", "e"]) => 90.0, Set(["c", "d"]) => 70.0)

        l = nj(d; ladderize=true)
        @test pathdist(l, names) ≈ P
    end
end

@testset "collections of name => sequence pairs" begin
    aln = readalignment(fixture("small.fa"))
    ps = collect(aln)                                     # name => sequence pairs
    gen = (n => s for (n, s) in ps)                       # any iterable of pairs
    @test Alignment(gen) == aln
    @test Alignment(aln) === aln
    @test distances(gen) == distances(aln)
    @test basefreqs(gen) == basefreqs(aln)
    @test writenewick(nj(gen; progress=false)) == writenewick(nj(aln; progress=false))
    @test writenewick(upgma(gen; progress=false)) == writenewick(upgma(aln; progress=false))
    @test_throws ArgumentError nj((1, 2, 3))
end
