def integer_cuberoot(n: int) -> int:
    if n < 0:
        raise ValueError("integer_cuberoot is only defined for nonnegative integers")

    lo, hi = 0, 1
    while hi ** 3 <= n:
        hi *= 2

    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid ** 3 <= n:
            lo = mid
        else:
            hi = mid
    return lo


def rho3(x: tuple[int, int, int]) -> int:
    a, b, c = x
    if a == b == c != 0:
        return (a + 1) ** 3 - 3
    elif b == c == 0:
        return (a + 1) ** 3 - 1
    elif a == b and c == 0:
        return (a + 1) ** 3 - 2
    elif a < c and b <= c:
        return c ** 3 + 3 * b * c + 3 * a
    elif b < a and c <= a:
        return a ** 3 + 3 * a * c + 3 * b - 1
    elif c < b and a <= b:
        return b ** 3 + 3 * a * b + 3 * ((c - 1) % b) + 1
    else:
        raise ValueError(f"Not a canonical order-3 onion state: {x}")


def anti_rho3(N: int) -> tuple[int, int, int]:
    if N < 0:
        raise ValueError("Ranks must be nonnegative")
    if N == 0:
        return 0, 0, 0

    m = integer_cuberoot(N)
    if N == (m + 1) ** 3 - 3:
        return m, m, m
    if N == (m + 1) ** 3 - 2:
        return m, m, 0
    if N == (m + 1) ** 3 - 1:
        return m, 0, 0

    t = N - m ** 3
    e = t % 3
    s = t // 3
    w = s % m
    v = s // m
    delta = (w + 1) % m

    if e == 0:
        return w, v, m
    if e == 1:
        return v, m, delta
    if delta == 0:
        return m, 0, v + 1
    return m, delta, v


def lambda3(x: tuple[int, int, int]) -> tuple[int, int, int]:
    N = rho3(x)
    m = integer_cuberoot(N)
    t = N - m ** 3
    return m, t // 3, t % 3


def omega3(m: int, u: int, e: int) -> tuple[int, int, int]:
    if m == 0:
        if (u, e) != (0, 0):
            raise ValueError(f"Invalid canonical order-3 data: {(m, u, e)}")
        return 0, 0, 0

    if e == 0:
        if u == m ** 2 + m:
            return m, 0, 0
        v, w = divmod(u, m)
        return w, v, m
    if e == 1:
        if u == m ** 2 + m - 1:
            return m, m, m
        v, w = divmod(u, m)
        return v, m, (w + 1) % m
    if e == 2:
        if u == m ** 2 + m - 1:
            return m, m, 0
        v, w = divmod(u, m)
        delta = (w + 1) % m
        if delta == 0:
            return m, 0, v + 1
        return m, delta, v
    raise ValueError(f"Invalid order-3 branch: {e}")


def carry_down(m: int, u: int, e: int) -> tuple[int, int, int]:
    while u < 0:
        if e == 0:
            m, u, e = m - 1, u + m ** 2 - m, 1
        elif e == 1:
            m, u, e = m - 1, u + m ** 2 - m, 2
        else:
            m, u, e = m - 1, u + m ** 2 - m + 1, 0
    return m, u, e


def add_states3(x: tuple[int, int, int], y: tuple[int, int, int]) -> tuple[int, int, int]:
    m, u, e = lambda3(x)
    n, v, t = lambda3(y)
    M = m + n
    U = u + v - m * n * (m + n) + (e + t) // 3
    E = (e + t) % 3
    M, U, E = carry_down(M, U, E)
    return omega3(M, U, E)


def carry_up(m: int, u: int, e: int) -> tuple[int, int, int]:
    while (e == 0 and u > m ** 2 + m) or (e in (1, 2) and u > m ** 2 + m - 1):
        if e == 0:
            m, u, e = m + 1, u - m ** 2 - m - 1, 2
        elif e == 1:
            m, u, e = m + 1, u - m ** 2 - m, 0
        else:
            m, u, e = m + 1, u - m ** 2 - m, 1
    return m, u, e


def mul_states3(x: tuple[int, int, int], y: tuple[int, int, int]) -> tuple[int, int, int]:
    m, u, e = lambda3(x)
    n, v, t = lambda3(y)
    M = m * n
    S = m ** 3 * t + n ** 3 * e + e * t
    U = m ** 3 * v + n ** 3 * u + 3 * u * v + u * t + v * e + S // 3
    E = S % 3
    M, U, E = carry_up(M, U, E)
    return omega3(M, U, E)


def divmod_states3(
    x: tuple[int, int, int], y: tuple[int, int, int]
) -> tuple[tuple[int, int, int], tuple[int, int, int]]:
    divisor = rho3(y)
    if divisor == 0:
        raise ZeroDivisionError("division by zero onion state")
    quotient, remainder = divmod(rho3(x), divisor)
    return anti_rho3(quotient), anti_rho3(remainder)


def floordiv_states3(x: tuple[int, int, int], y: tuple[int, int, int]) -> tuple[int, int, int]:
    return divmod_states3(x, y)[0]


def mod_states3(x: tuple[int, int, int], y: tuple[int, int, int]) -> tuple[int, int, int]:
    return divmod_states3(x, y)[1]
