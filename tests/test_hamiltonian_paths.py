from itertools import product
from math import factorial, isqrt, prod
import pytest


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


@pytest.mark.parametrize(('order', 'top_layer'), ((2, 5), (3, 3), (4, 2)))
def test_small_universal_properties(order: int, top_layer: int) -> int:
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


@pytest.mark.parametrize(('n', 'k'), tuple(product(range(2, 6), range(2, 6))))
def test_onion_theorem(n: int, k: int):
    if k ** n > 625:
        return
    current = prefer_max_words(k, n)
    previous = prefer_max_words(k - 1, n)
    assert current[-(k - 1) ** n :] == previous
    reversed_current = reverse_prefer_max_words(k, n)
    reversed_previous = reverse_prefer_max_words(k - 1, n)
    assert reversed_current[: (k - 1) ** n] == reversed_previous
    assert_debruijn(reversed_current, k, n)


@pytest.mark.parametrize(('n', 'k'), ((1, 1), (1, 2), (1, 5), (2, 2), (2, 3), (2, 4), (2, 5), (3, 2), (3, 3), (4, 2)))
def test_enumeration_layers(n: int, k: int):
    brute_force = count_layer_hamiltonian_cycles(n, k)
    closed_form = layer_cycle_formula(n, k)
    assert brute_force == closed_form, ((n, k), brute_force, closed_form)


@pytest.mark.parametrize(('n', 'j'), ((1, 1), (1, 5), (2, 2), (2, 3), (2, 4), (2, 5), (3, 2), (3, 3), (4, 2)))
def test_enumeration_prefixes(n: int, j: int):
    brute_force = sum(1 for _ in iter_finite_onion_prefixes(n, j))
    closed_form = onion_prefix_formula(n, j)
    assert brute_force == closed_form, ((n, j), brute_force, closed_form)
