from itertools import product
from math import factorial, isqrt, prod

from order2_debruijn_arithmetic import (
    add_states,
    anti_rho2,
    lambda2,
    mul_states,
    no_sqrt_rho2,
    rho2,
)
from order3_debruijn_arithmetic import (
    add_states3,
    anti_rho3,
    integer_cuberoot,
    lambda3,
    mul_states3,
    omega3,
    rho3,
)


def prefer_max_words(k: int, n: int) -> list[tuple[int, ...]]:
    word = (0,) * (n - 1) + (k - 1,)
    words = [word]
    seen = {word}

    while len(words) < k ** n:
        prefix = word[1:]
        for tau in range(k - 1, -1, -1):
            candidate = prefix + (tau,)
            if candidate not in seen:
                words.append(candidate)
                seen.add(candidate)
                word = candidate
                break
        else:
            raise RuntimeError(f"prefer-max got stuck for {(k, n)}")
    return words


def reverse_prefer_max_words(k: int, n: int) -> list[tuple[int, ...]]:
    return [tuple(reversed(word)) for word in reversed(prefer_max_words(k, n))]


def onion_prefix(order: int, top_layer: int) -> list[tuple[int, ...]]:
    return reverse_prefer_max_words(top_layer, order)


def assert_debruijn(words: list[tuple[int, ...]], k: int, n: int) -> None:
    assert len(words) == k ** n, f"Wrong length for {(k, n)}"
    assert set(words) == set(product(range(k), repeat=n)), f"Wrong vertex set for {(k, n)}"
    for left, right in zip(words, words[1:] + words[:1]):
        assert left[1:] == right[:-1], f"Broken overlap in {(k, n)}: {left} -> {right}"


def choose_layer_for_rank(order: int, rank: int) -> int:
    k = 1
    while k ** order <= rank:
        k += 1
    return k


def count_order2_add_carries(x: tuple[int, int], y: tuple[int, int]) -> int:
    m, u, e = lambda2(x)
    n, v, t = lambda2(y)
    m, u, e = m + n, u + v - m * n - (e + t) // 2, (e + t) % 2
    carries = 0
    while (e == 0 and u < 0) or (e == 1 and u <= 0):
        carries += 1
        if e == 0:
            m, u, e = m - 1, u + m, 1
        else:
            m, u, e = m - 1, u + m - 1, 0
    return carries


def count_order2_mul_carries(x: tuple[int, int], y: tuple[int, int]) -> int:
    m, u, e = lambda2(x)
    n, v, t = lambda2(y)
    s = t * m ** 2 + e * n ** 2 - e * t
    m, u, e = (
        m * n,
        m ** 2 * v + n ** 2 * u + 2 * u * v - u * t - v * e - s // 2,
        s % 2,
    )
    carries = 0
    while u > m:
        carries += 1
        m, u, e = m + 1, u - m - e, 1 - e
    return carries


def count_order3_add_carries(x: tuple[int, int, int], y: tuple[int, int, int]) -> int:
    m, u, e = lambda3(x)
    n, v, t = lambda3(y)
    m, u, e = m + n, u + v - m * n * (m + n) + (e + t) // 3, (e + t) % 3
    carries = 0
    while u < 0:
        carries += 1
        if e == 0:
            m, u, e = m - 1, u + m ** 2 - m, 1
        elif e == 1:
            m, u, e = m - 1, u + m ** 2 - m, 2
        else:
            m, u, e = m - 1, u + m ** 2 - m + 1, 0
    return carries


def count_order3_mul_carries(x: tuple[int, int, int], y: tuple[int, int, int]) -> int:
    m, u, e = lambda3(x)
    n, v, t = lambda3(y)
    s = m ** 3 * t + n ** 3 * e + e * t
    m, u, e = (
        m * n,
        m ** 3 * v + n ** 3 * u + 3 * u * v + u * t + v * e + s // 3,
        s % 3,
    )
    carries = 0
    while (e == 0 and u > m ** 2 + m) or (e in (1, 2) and u > m ** 2 + m - 1):
        carries += 1
        if e == 0:
            m, u, e = m + 1, u - m ** 2 - m - 1, 2
        elif e == 1:
            m, u, e = m + 1, u - m ** 2 - m, 0
        else:
            m, u, e = m + 1, u - m ** 2 - m, 1
    return carries


def reverse_prefer_max_prefix_by_rank(order: int, max_rank: int) -> list[tuple[int, ...]]:
    return reverse_prefer_max_words(choose_layer_for_rank(order, max_rank), order)


def layer_vertices(n: int, k: int) -> list[tuple[int, ...]]:
    return [word for word in product(range(k), repeat=n) if (k - 1) in word]


def layer_adjacency(n: int, k: int) -> dict[tuple[int, ...], list[tuple[int, ...]]]:
    adjacency: dict[tuple[int, ...], list[tuple[int, ...]]] = {}
    for word in layer_vertices(n, k):
        suffix = word[1:]
        adjacency[word] = []
        for tau in range(k):
            candidate = suffix + (tau,)
            if (k - 1) in candidate:
                adjacency[word].append(candidate)
    return adjacency


def count_layer_hamiltonian_cycles(n: int, k: int) -> int:
    if n == 1:
        return 1

    adjacency = layer_adjacency(n, k)
    start = min(adjacency)
    seen = {start}
    vertex_count = len(adjacency)

    def dfs(vertex: tuple[int, ...]) -> int:
        if len(seen) == vertex_count:
            return int(start in adjacency[vertex])

        total = 0
        for nxt in adjacency[vertex]:
            if nxt not in seen:
                seen.add(nxt)
                total += dfs(nxt)
                seen.remove(nxt)
        return total

    return dfs(start)


def list_layer_paths(n: int, k: int) -> list[list[tuple[int, ...]]]:
    if n == 1:
        return [[(k - 1,)]]

    adjacency = layer_adjacency(n, k)
    start = (0,) * (n - 1) + (k - 1,)
    end = (k - 1,) + (0,) * (n - 1)
    seen = {start}
    path = [start]
    paths: list[list[tuple[int, ...]]] = []
    vertex_count = len(adjacency)

    def dfs(vertex: tuple[int, ...]) -> None:
        if len(path) == vertex_count:
            if vertex == end:
                paths.append(path.copy())
            return

        for nxt in adjacency[vertex]:
            if nxt not in seen:
                seen.add(nxt)
                path.append(nxt)
                dfs(nxt)
                path.pop()
                seen.remove(nxt)

    dfs(start)
    return paths


def layer_cycle_formula(n: int, k: int) -> int:
    exponent = k ** (n - 1) - (k - 1) ** (n - 1)
    return factorial(k) ** exponent // (k ** (n - 1))


def onion_prefix_formula(n: int, j: int) -> int:
    return prod(layer_cycle_formula(n, k) for k in range(2, j + 1))


def iter_finite_onion_prefixes(n: int, j: int):
    base = [(0,) * n]
    if j == 1:
        yield base
        return

    layer_paths = {k: list_layer_paths(n, k) for k in range(2, j + 1)}

    def extend(k: int, prefix: list[tuple[int, ...]]):
        if k > j:
            yield prefix
            return
        for path in layer_paths[k]:
            yield from extend(k + 1, prefix + path)

    yield from extend(2, base)


def verify_small_universal_properties(order: int, top_layer: int) -> int:
    checked = 0
    for prefix in iter_finite_onion_prefixes(order, top_layer):
        checked += 1
        for k in range(1, top_layer + 1):
            assert_debruijn(prefix[: k ** order], k, order)

        rank = {word: index for index, word in enumerate(prefix)}

        for index, word in enumerate(prefix):
            mu = max(word)
            assert mu ** order <= index < (mu + 1) ** order
            for step in range(len(prefix) - index):
                new_word = prefix[index + step]
                t = max(new_word)
                assert t ** order <= index + step < (t + 1) ** order

        if order == 1:
            continue

        alphabet = range(top_layer)
        for context in product(alphabet, repeat=order - 1):
            start = max(context)
            indices = [rank[context + (tau,)] for tau in range(start, top_layer)]
            assert indices == sorted(indices)

        for left_len in range(order):
            right_len = order - 1 - left_len
            for left in product(alphabet, repeat=left_len):
                for right in product(alphabet, repeat=right_len):
                    context = left + right
                    mu = max(context) if context else -1
                    indices = [rank[left + (sigma,) + right] for sigma in range(mu + 1, top_layer)]
                    assert indices == sorted(indices)
                    for sigma in range(mu + 1, top_layer):
                        word = left + (sigma,) + right
                        assert sigma ** order <= rank[word] < (sigma + 1) ** order
    return checked


def verify_onion_theorem() -> int:
    checked = 0
    for n in range(2, 6):
        for k in range(2, 6):
            if k ** n > 625:
                continue
            current = prefer_max_words(k, n)
            previous = prefer_max_words(k - 1, n)
            assert current[-(k - 1) ** n :] == previous
            reversed_current = reverse_prefer_max_words(k, n)
            reversed_previous = reverse_prefer_max_words(k - 1, n)
            assert reversed_current[: (k - 1) ** n] == reversed_previous
            assert_debruijn(reversed_current, k, n)
            checked += 1
    return checked


def verify_enumeration_results() -> tuple[int, int]:
    layer_cases = [(1, 1), (1, 2), (1, 5), (2, 2), (2, 3), (2, 4), (2, 5), (3, 2), (3, 3), (4, 2)]
    for n, k in layer_cases:
        brute_force = count_layer_hamiltonian_cycles(n, k)
        closed_form = layer_cycle_formula(n, k)
        assert brute_force == closed_form, ((n, k), brute_force, closed_form)

    prefix_cases = [(1, 1), (1, 5), (2, 2), (2, 3), (2, 4), (2, 5), (3, 2), (3, 3), (4, 2)]
    for n, j in prefix_cases:
        brute_force = sum(1 for _ in iter_finite_onion_prefixes(n, j))
        closed_form = onion_prefix_formula(n, j)
        assert brute_force == closed_form, ((n, j), brute_force, closed_form)

    return len(layer_cases), len(prefix_cases)


def verify_order2_results() -> tuple[int, int]:
    words = reverse_prefer_max_prefix_by_rank(2, 600)
    rank = {word: index for index, word in enumerate(words)}

    for index, word in enumerate(words[:500]):
        assert rho2(word) == index
        assert anti_rho2(index) == word
        assert no_sqrt_rho2(word) == index

    for k in range(1, 16):
        layer = words[k ** 2 : (k + 1) ** 2]
        expected = [(0, k)]
        for j in range(1, k):
            expected.extend([(k, j), (j, k)])
        expected.extend([(k, k), (k, 0)])
        assert layer == expected, (k, layer, expected)
        for j in range(k):
            assert rho2((j, k)) == k ** 2 + 2 * j
        for j in range(1, k + 1):
            assert rho2((k, j)) == k ** 2 + 2 * j - 1
        assert rho2((k, 0)) == (k + 1) ** 2 - 1

    for left in words[:120]:
        for right in words[:120]:
            left_layer = lambda2(left)[0]
            right_layer = lambda2(right)[0]
            add_rank = rho2(left) + rho2(right)
            mul_rank = rho2(left) * rho2(right)
            add_carries = count_order2_add_carries(left, right)
            mul_carries = count_order2_mul_carries(left, right)

            assert add_states(left, right) == anti_rho2(add_rank)
            assert mul_states(left, right) == anti_rho2(mul_rank)
            assert add_carries == left_layer + right_layer - isqrt(add_rank)
            assert mul_carries == isqrt(mul_rank) - left_layer * right_layer
            assert add_carries <= left_layer + right_layer
            assert mul_carries <= left_layer + right_layer

    assert words[:16] == [
        (0, 0),
        (0, 1),
        (1, 1),
        (1, 0),
        (0, 2),
        (2, 1),
        (1, 2),
        (2, 2),
        (2, 0),
        (0, 3),
        (3, 1),
        (1, 3),
        (3, 2),
        (2, 3),
        (3, 3),
        (3, 0),
    ]
    assert [rank[word] for word in [(1, 2), (2, 0), (0, 3), (3, 1)]] == [6, 8, 9, 10]
    assert [rank[word] for word in [(1, 0), (2, 0), (3, 0), (4, 0)]] == [3, 8, 15, 24]
    assert [rank[word] for word in [(1, 2), (1, 3), (1, 4)]] == [6, 11, 18]
    assert add_states((2, 4), (3, 1)) == (5, 3)
    assert mul_states((1, 2), (2, 0)) == (6, 0)

    return 500, 120 ** 2


def verify_order3_results() -> tuple[int, int]:
    words = reverse_prefer_max_prefix_by_rank(3, 1200)
    rank = {word: index for index, word in enumerate(words)}

    for index, word in enumerate(words[:900]):
        assert rho3(word) == index
        assert anti_rho3(index) == word
        assert omega3(*lambda3(word)) == word

    for left in words[:90]:
        for right in words[:90]:
            left_layer = lambda3(left)[0]
            right_layer = lambda3(right)[0]
            add_rank = rho3(left) + rho3(right)
            mul_rank = rho3(left) * rho3(right)
            add_carries = count_order3_add_carries(left, right)
            mul_carries = count_order3_mul_carries(left, right)

            assert add_states3(left, right) == anti_rho3(add_rank)
            assert mul_states3(left, right) == anti_rho3(mul_rank)
            assert add_carries == left_layer + right_layer - integer_cuberoot(add_rank)
            assert mul_carries == integer_cuberoot(mul_rank) - left_layer * right_layer
            assert add_carries <= left_layer + right_layer
            assert mul_carries <= left_layer + right_layer

    assert words[8:27] == [
        (0, 0, 2),
        (0, 2, 1),
        (2, 1, 0),
        (1, 0, 2),
        (0, 2, 0),
        (2, 0, 1),
        (0, 1, 2),
        (1, 2, 1),
        (2, 1, 1),
        (1, 1, 2),
        (1, 2, 0),
        (2, 0, 2),
        (0, 2, 2),
        (2, 2, 1),
        (2, 1, 2),
        (1, 2, 2),
        (2, 2, 2),
        (2, 2, 0),
        (2, 0, 0),
    ]
    assert rho3((1, 2, 0)) == 18
    assert anti_rho3(23) == (1, 2, 2)
    assert add_states3((1, 2, 0), (0, 2, 1)) == (0, 0, 3)

    return 900, 90 ** 2


def main() -> None:
    onion_cases = verify_onion_theorem()
    layer_cases, prefix_cases = verify_enumeration_results()
    small_prefixes = 0
    for order, top_layer in [(2, 5), (3, 3), (4, 2)]:
        small_prefixes += verify_small_universal_properties(order, top_layer)
    order2_states, order2_pairs = verify_order2_results()
    order3_states, order3_pairs = verify_order3_results()

    print("Verified onion theorem on", onion_cases, "small (n, k) instances.")
    print("Verified layer-count formula on", layer_cases, "small cases.")
    print("Verified compatible onion-prefix count on", prefix_cases, "small cases.")
    print("Verified universal structural/arithmetic consequences on", small_prefixes, "exhaustively generated finite onion prefixes.")
    print("Verified order-2 rank/unrank on", order2_states, "states and direct addition/multiplication plus carry counts on", order2_pairs, "pairs.")
    print("Verified order-3 rank/unrank on", order3_states, "states and direct addition/multiplication plus carry counts on", order3_pairs, "pairs.")


if __name__ == "__main__":
    main()
