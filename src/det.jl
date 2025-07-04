"""
    det_4x4(M::AbstractMatrix)

Compute the determinant of a 4x4 matrix using an explicit cofactor
expansion.
"""
function det_4x4(M::AbstractMatrix)
    size(M) == (4, 4) || throw(ArgumentError("M must be a 4x4 matrix"))

    # IMPROVE: Use muladd

    fmms(a, b, c, d) = a * b - c * d

    res =
        fmms(M[1, 4], M[2, 3], M[1, 3], M[2, 4]) * fmms(M[3, 2], M[4, 1], M[3, 1], M[4, 2])

    res +=
        fmms(M[1, 2], M[2, 4], M[1, 4], M[2, 2]) * fmms(M[3, 3], M[4, 1], M[3, 1], M[4, 3])

    res +=
        fmms(M[1, 3], M[2, 2], M[1, 2], M[2, 3]) * fmms(M[3, 4], M[4, 1], M[3, 1], M[4, 4])

    res +=
        fmms(M[1, 4], M[2, 1], M[1, 1], M[2, 4]) * fmms(M[3, 3], M[4, 2], M[3, 2], M[4, 3])

    res +=
        fmms(M[1, 1], M[2, 3], M[1, 3], M[2, 1]) * fmms(M[3, 4], M[4, 2], M[3, 2], M[4, 4])

    res +=
        fmms(M[1, 2], M[2, 1], M[1, 1], M[2, 2]) * fmms(M[3, 4], M[4, 3], M[3, 3], M[4, 4])

    return res
end
