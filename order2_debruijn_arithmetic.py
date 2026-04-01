def rho2(x: tuple[int, int]) -> int:
    a, b = x
    if a < b:
        return b ** 2 + 2 * a
    elif 0 < b <= a:
        return a ** 2 + 2 * b - 1
    else:
        return a ** 2 + 2 * a


def anti_rho2(N: int) -> tuple[int, int]:
    m = int(N ** 0.5)
    while (m + 1) ** 2 <= N:
        m += 1
    while m ** 2 > N:
        m -= 1
    t = N - m ** 2
    if t % 2:
        return m, (t + 1) // 2
    elif t < 2 * m:
        return t // 2, m
    else:
        return m, 0


def lambda2(x: tuple[int, int]) -> tuple[int, int, int]:
    a, b = x
    if a < b:
        return b, a, 0
    elif b == 0:
        return a, a, 0
    else:
        return a, b, 1


def no_sqrt_rho2(x: tuple[int, int]) -> int:
    m, u, e = lambda2(x)
    return m ** 2 + 2 * u - e


def omega2(m: int, u: int, e: int) -> tuple[int, int]:
    if e == 0:
        if u < m:
            return u, m
        return m, 0
    return m, u


def carry_down(m: int, u: int, e: int) -> tuple[int, int, int]:
    while (e == 0 and u < 0) or (e == 1 and u <= 0):
        if e == 0:
            m, u, e = m - 1, u + m, 1
        else:
            m, u, e = m - 1, u + m - 1, 0
    return m, u, e


def add_states(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    m, u, e = lambda2(x)
    n, v, t = lambda2(y)
    M = m + n
    U = u + v - m * n - (e + t) // 2
    E = (e + t) % 2
    M, U, E = carry_down(M, U, E)
    return omega2(M, U, E)


def carry_up(m: int, u: int, e: int) -> tuple[int, int, int]:
    while u > m:
        m, u, e = m + 1, u - m - e, 1 - e
    return m, u, e


def mul_states(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    m, u, e = lambda2(x)
    n, v, t = lambda2(y)
    M = m * n
    S = t * m ** 2 + e * n ** 2 - e * t
    U = m ** 2 * v + n ** 2 * u + 2 * u * v - u * t - v * e - S // 2
    E = S % 2
    M, U, E = carry_up(M, U, E)
    return omega2(M, U, E)


def divmod_states(x: tuple[int, int], y: tuple[int, int]) -> tuple[tuple[int, int], tuple[int, int]]:
    divisor = rho2(y)
    if divisor == 0:
        raise ZeroDivisionError("division by zero onion state")
    quotient, remainder = divmod(rho2(x), divisor)
    return anti_rho2(quotient), anti_rho2(remainder)


def floordiv_states(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    return divmod_states(x, y)[0]


def mod_states(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    return divmod_states(x, y)[1]
