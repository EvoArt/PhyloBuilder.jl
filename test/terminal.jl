@testset "progress bars and terminal styling" begin
    P = PhyloBuilder
    P.setprogress!(true); @test P._progress_default()
    P.setprogress!(false); @test !P._progress_default()
    P._PROGRESS[] = nothing
    withenv("CI" => "true") do
        @test !P._progress_default()
    end
    withenv("CI" => nothing) do
        @test P._progress_default() isa Bool
    end
    for n in (5, 0), color in (false, true)
        buf = IOBuffer()
        io = IOContext(buf, :color => color, :displaysize => (24, 100))
        p = P._progress(n, "Testing"; unit="items", enabled=true, io=io)
        p.tty = true                                   # draw as on a terminal
        P._render(p)
        p.tlast = 0.0; P._next!(p; info="first")
        p.tlast = 0.0; P._next!(p, 2)
        P._finish!(p; info="done", ok=color)
        P._finish!(p)                                  # a second finish does nothing
        out = String(take!(buf))
        @test occursin("Testing", out) && occursin("done", out) && occursin("\r", out)
        @test occursin("\e[38;5;", out) == color
    end
    buf = IOBuffer()
    p = P._progress(3, "Quiet"; enabled=true, io=buf)
    P._next!(p, 2); P._finish!(p; ok=false)
    @test count(==('\n'), String(take!(buf))) == 1     # not a terminal: the final line only
    q = P._progress(3, "Off"; enabled=false, io=devnull)
    @test P._next!(q) === q && P._finish!(q) === q
    @test P._secs(0.5) == "500 ms" && P._secs(5.0) == "5.0 s"
    @test P._secs(125.0) == "2 min 5 s" && P._secs(7300.0) == "2 h 2 min"
    @test P._fmt(-1234567) == "-1,234,567" && P._fmt(3.0) == "3" && P._fmt(0.25) == "0.25"
    @test P._g(NaN) == "NaN" && P._g(Inf) == "Inf" && P._g(-Inf) == "-Inf" && P._g(0.0) == "0"
    @test P._g(12345.0) == "1.23e+04"
    @test P._sty(IOContext(IOBuffer(), :color => true), "x"; bold=true) == "\e[1mx\e[0m"
    @test P._sty(IOContext(IOBuffer(), :color => true), "x") == "x"
    # a real distance matrix of 200+ taxa draws its bar, saturated pairs noted
    rng = MersenneTwister(5)
    big = Alignment(["t$i" for i in 1:200], [join(rand(rng, "ACGT", 40)) for _ in 1:200])
    err = mktemp() do path, f
        redirect_stderr(f) do
            distances(big; model=:JC69, progress=true)
        end
        close(f); read(path, String)
    end
    @test occursin("Distances", err) && occursin("saturated", err)
end
