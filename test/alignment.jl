@testset "Alignment" begin
    @testset "construction" begin
        a = Alignment(["x", "y", "z"], ["ACGT-", "acgtN", "RYAC."])
        @test ntaxa(a) == 3 && nsites(a) == 5 && size(a) == (3, 5) && length(a) == 3
        @test taxa(a) == ["x", "y", "z"]
        @test a[1] == "ACGT-"
        @test a["y"] == "ACGTN"              # upper case
        @test a[3] == "RYAC-"                # '.' is a gap
        @test sequences(a) == ["ACGT-", "ACGTN", "RYAC-"]
        @test collect(a)[2] == ("y" => "ACGTN")
        @test Alignment(["x" => "AC", "y" => "GT"]) == Alignment(["x", "y"], ["AC", "GT"])
        @test taxa(Alignment(["AC", "GT"])) == ["seq1", "seq2"]
        @test sort(taxa(Alignment(Dict("p" => "AC", "q" => "GT")))) == ["p", "q"]
        @test Alignment(["x"], ["AC GT\n"])[1] == "ACGT"   # whitespace ignored
        @test Alignment(["x"], ["ACGU"])[1] == "ACGT"
        @test_throws ArgumentError Alignment(["x", "y"], ["ACG", "AC"])
        @test_throws ArgumentError Alignment(["x", "x"], ["AC", "AC"])
        @test_throws ArgumentError Alignment(["x"], ["ACZ"])
        @test_throws ArgumentError Alignment(["x", "y"], ["AC"])
        @test_throws KeyError a["nope"]
    end

    @testset "subsetting" begin
        a = Alignment(["x", "y", "z"], ["ACGT", "CCGG", "TTAA"])
        @test a[[1, 3], 2:3] == Alignment(["x", "z"], ["CG", "TA"])
        @test a[["z", "x"], :] == Alignment(["z", "x"], ["TTAA", "ACGT"])
        @test a[[true, false, true], 4] == Alignment(["x", "z"], ["T", "A"])
        @test a["y", :] == Alignment(["y"], ["CCGG"])
        @test copy(a) == a
    end

    @testset "base frequencies" begin
        a = Alignment(["x", "y"], ["AAR-", "CGTN"])
        # R = A|G shared out; gap and N ignored
        @test basefreqs(a) ≈ [2.5, 1, 1.5, 1] ./ 6
        @test basefreqs(a; ambiguous=false) ≈ [2, 1, 1, 1] ./ 5
        @test basefreqs(Alignment(["x"], ["NN--"])) == fill(0.25, 4)
        # ape::base.freq counts only A, C, G and T
        ref = parse.(Float64, readlines(fixture("basefreq.txt")))
        @test basefreqs(readalignment(fixture("small.fa")); ambiguous=false) ≈ ref atol = 1e-15
    end

    @testset "site patterns" begin
        a = Alignment(["x", "y", "z"], ["AACA-", "AACAN", "GGTGA"])
        p = PhyloBuilder._patterns(a)
        @test PhyloBuilder._npatterns(p) == 3      # AAG ×3, CCT, (gap)(N)A — gap and N both missing
        @test sum(p.weight) == nsites(a)
        @test p.data[:, p.site] == a.data .& 0x0f
        w = PhyloBuilder._resample(p, Xoshiro(1))
        @test sum(w) == nsites(a) && all(>=(0), w)
        @test PhyloBuilder._siteclass(UInt8[0x01, 0x01, 0x05]) == (true, false)   # A A R: constant
        @test PhyloBuilder._siteclass(UInt8[0x01, 0x01, 0x02, 0x02]) == (false, true)
        @test PhyloBuilder._siteclass(UInt8[0x01, 0x01, 0x02, 0x0f]) == (false, false)
    end

    @testset "reading and writing" begin
        a = Alignment(["s1", "long_name_2", "s3"], ["ACGTRYN-ACGT", "ACGTACGTACGT", "TTTT----CCCC"])
        mktempdir() do dir
            f = writefasta(joinpath(dir, "a.fa"), a; width=5)
            @test readalignment(f) == a
            @test readfasta(f) == a
            @test count(==('\n'), read(f, String)) == 3 * (1 + 3)
            f = writephylip(joinpath(dir, "a.phy"), a)
            @test readalignment(f) == a
            @test readphylip(f) == a

            f = joinpath(dir, "inter.phy")     # interleaved PHYLIP
            write(f, "3 12\ns1  ACGTRY\nlong_name_2 ACGTAC\ns3  TTTT--\nN-ACGT\nGTACGT\n--CCCC\n")
            @test readalignment(f) == a

            f = joinpath(dir, "desc.fa")       # FASTA descriptions dropped when ids are unique
            write(f, ">a some description\nAC\nGT\n>b other\nACGG\n")
            @test taxa(readalignment(f)) == ["a", "b"]
            write(f, ">a one\nACGT\n>a two\nACGG\n")
            @test taxa(readalignment(f)) == ["a one", "a two"]

            f = joinpath(dir, "a.nex")
            write(f, """
                #NEXUS
                begin data;
                  dimensions ntax=2 nchar=8;
                  format datatype=dna missing=? gap=-;
                  matrix
                    'taxon one' ACGT [comment] ACG-
                    t2          ACGTACGT
                    'taxon one'
                  ;
                end;
                """)
            n = readalignment(f)
            @test taxa(n) == ["taxon one", "t2"]
            @test sequences(n) == ["ACGTACG-", "ACGTACGT"]

            f = joinpath(dir, "bad.txt")
            write(f, "hello")
            @test_throws ArgumentError readalignment(f)
        end
    end
end
