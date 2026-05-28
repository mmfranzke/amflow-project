import mpmath as mp


# Small numerical utilities used by scripts and tests.


def relative_error(value, reference):
    # Use an absolute error if the reference value is exactly zero.
    if abs(reference) == 0:
        return abs(value - reference)
    return abs((value - reference) / reference)


def quad_split(function, a, b, parts=4):
    # Split intervals to improve stability near endpoint singularities.
    points = [a + (b - a) * i / parts for i in range(parts + 1)]
    return mp.quad(function, points)

