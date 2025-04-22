@testset "verify_and_refine_root" begin
    f = x -> [sin(x[1]), sin(x[2]), sin(x[3]), cos(x[4])]
    df = x -> [cos(x[1]) 0 0 0; 0 cos(x[2]) 0 0; 0 0 cos(x[3]) 0; 0 0 0 -sin(x[4])]

    root_true = SVector(Arb(0), Arb(π), 2Arb(π), Arb(π) / 2)
    root_enclosure = add_error.(root_true, Mag(1))

    root = CGL2.verify_and_refine_root(f, df, root_enclosure)

    atol = 0
    rtol = 4eps(one(first(root)))

    @test all(Arblib.overlaps.(root, root_true))
    @test all(ArbExtras.check_tolerance.(root; atol, rtol))
end

@testset "verify_root_from_approximation" begin
    f = x -> [sin(x[1]), sin(x[2]), sin(x[3]), cos(x[4])]
    df = x -> [cos(x[1]) 0 0 0; 0 cos(x[2]) 0 0; 0 0 cos(x[3]) 0; 0 0 0 -sin(x[4])]

    root_true = SVector(Arb(0), Arb(π), 2Arb(π), Arb(π) / 2)
    root_approximation = root_true .+ 0.1

    root, root_uniqueness = CGL2.verify_root_from_approximation(f, df, root_approximation)

    atol = 0
    rtol = 4eps(one(first(root)))

    @test all(Arblib.overlaps.(root, root_true))
    @test all(ArbExtras.check_tolerance.(root; atol, rtol))
    @test all(Arblib.contains_interior.(root_uniqueness, root))

    complex_f = x -> [exp(x[1]) + exp(Arb(1)), sin(x[2]), sin(x[3]), cos(x[4])]
    complex_df = x -> [exp(x[1]) 0 0 0; 0 cos(x[2]) 0 0; 0 0 cos(x[3]) 0; 0 0 0 -sin(x[4])]

    complex_root_true = SVector(Acb(1, π), Acb(π), 2Acb(π), Acb(π) / 2)
    complex_root_approximation = complex_root_true .+ (0.1 + 0.2im)

    complex_root, complex_root_uniqueness = CGL2.verify_root_from_approximation(
        complex_f,
        complex_df,
        complex_root_approximation,
    )

    atol = eps(one(first(root)))
    rtol = 8eps(one(first(root)))

    @test all(Arblib.overlaps.(complex_root, complex_root_true))
    @test all(ArbExtras.check_tolerance.(real.(complex_root); atol, rtol))
    @test all(ArbExtras.check_tolerance.(imag.(complex_root); atol, rtol))
    @test all(Arblib.contains_interior.(complex_root_uniqueness, complex_root))
end
