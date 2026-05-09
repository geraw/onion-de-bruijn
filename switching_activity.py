from collections import Counter
from math import ceil, log2


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


def tail(word: tuple[int, ...]) -> bool:
    if word[-1] == 0:
        return False

    ell = 0
    while ell < len(word) and word[ell] == 0:
        ell += 1
    if ell == len(word):
        return False

    canonical = word[ell:] + word[:ell]
    return canonical == max(word[i:] + word[:i] for i in range(len(word)))


def reverse_tail_candidates(x: tuple[int, ...], k: int) -> list[int]:
    reversed_x = tuple(reversed(x))
    return [tau for tau in range(1, k) if tail(reversed_x + (tau,))]


def onion_successor_rule(word: tuple[int, ...]) -> tuple[int, ...]:
    sigma = word[0]
    x = word[1:]
    k = max(word) + 2
    reversed_x = tuple(reversed(x))
    candidates = reverse_tail_candidates(x, k)

    append_candidates: list[int] = []

    # Inverse of the forward prefer-max tail descent branch.
    if sigma + 1 < k and tail(reversed_x + (sigma + 1,)):
        append_candidates.append(sigma + 1)

    # Inverse of the forward prefer-max branch that appends max(T).
    if candidates and max(candidates) == sigma:
        append_candidates.append(0)

    # Inverse of the forward pass-through branch.
    if sigma == 0:
        if not candidates:
            append_candidates.append(0)
    else:
        if not tail(reversed_x + (sigma,)):
            append_candidates.append(sigma)

    if len(append_candidates) != 1:
        raise AssertionError(f"Non-unique onion successor for {word}: {append_candidates}")

    return x + (append_candidates[0],)


def bits(value: int, width: int) -> tuple[int, ...]:
    return tuple((value >> i) & 1 for i in range(width))


def gray(value: int) -> int:
    return value ^ (value >> 1)


def hamming(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return sum(a != b for a, b in zip(left, right))


def binary_counter_switching_stats(num_states: int) -> dict[str, object]:
    width = ceil(log2(num_states))
    bit_histogram: Counter[int] = Counter()

    for state in range(num_states):
        nxt = 0 if state == num_states - 1 else state + 1
        bit_histogram[hamming(bits(state, width), bits(nxt, width))] += 1

    average_toggles = sum(toggles * count for toggles, count in bit_histogram.items()) / num_states

    return {
        "num_states": num_states,
        "bit_width": width,
        "bit_histogram": dict(sorted(bit_histogram.items())),
        "average_bit_toggles": average_toggles,
        "max_bit_toggles": max(bit_histogram),
    }


def onion_words_by_successor(order: int, max_symbol: int) -> list[tuple[int, ...]]:
    reference_words = reverse_prefer_max_words(max_symbol + 1, order)
    num_states = len(reference_words)
    word = (0,) * order
    assert word == reference_words[0]
    words = [word]
    seen = {word}

    for step, expected in enumerate(reference_words[1:], start=1):
        word = onion_successor_rule(word)
        assert max(word) <= max_symbol, f"Successor escaped bounded layer set: {word}"
        assert word not in seen, f"Successor repeated too early: {word}"
        assert word == expected, f"Successor diverged from bounded reverse prefer-max orbit at step {step}: expected {expected}, got {word}"
        seen.add(word)
        words.append(word)

    # Iterating the current onion successor from 0^n must reproduce the bounded
    # reverse prefer-max orbit on [max_symbol + 1]^order exactly.
    assert words == reference_words
    return words


def onion_pointer_switching_stats(order: int, max_symbol: int, gray_pointer: bool) -> dict[str, object]:
    words = onion_words_by_successor(order, max_symbol)
    num_states = len(words)
    symbol_width = max(1, ceil(log2(max_symbol + 1)))
    pointer_width = max(1, ceil(log2(order)))

    buffer = list(words[0])
    head = 0

    symbol_bit_histogram: Counter[int] = Counter()
    pointer_bit_histogram: Counter[int] = Counter()
    total_bit_histogram: Counter[int] = Counter()
    changed_symbol_cell_histogram: Counter[int] = Counter()

    for index, word in enumerate(words):
        nxt = words[(index + 1) % num_states]
        assert nxt[:-1] == word[1:]

        new_symbol = nxt[-1]
        old_symbol = buffer[head]
        symbol_toggles = hamming(bits(old_symbol, symbol_width), bits(new_symbol, symbol_width))

        next_head = (head + 1) % order
        current_pointer = gray(head) if gray_pointer else head
        next_pointer = gray(next_head) if gray_pointer else next_head
        pointer_toggles = hamming(bits(current_pointer, pointer_width), bits(next_pointer, pointer_width))

        symbol_bit_histogram[symbol_toggles] += 1
        pointer_bit_histogram[pointer_toggles] += 1
        total_bit_histogram[symbol_toggles + pointer_toggles] += 1
        changed_symbol_cell_histogram[int(old_symbol != new_symbol)] += 1

        # Physical circular-buffer update: overwrite one location and advance the head.
        buffer[head] = new_symbol
        head = next_head

        logical_word = tuple(buffer[(head + offset) % order] for offset in range(order))
        assert logical_word == nxt

    average_symbol_toggles = sum(toggles * count for toggles, count in symbol_bit_histogram.items()) / num_states
    average_pointer_toggles = sum(toggles * count for toggles, count in pointer_bit_histogram.items()) / num_states
    average_total_toggles = sum(toggles * count for toggles, count in total_bit_histogram.items()) / num_states

    return {
        "num_states": num_states,
        "order": order,
        "max_symbol": max_symbol,
        "symbol_width": symbol_width,
        "pointer_width": pointer_width,
        "pointer_encoding": "gray" if gray_pointer else "binary",
        "symbol_bit_histogram": dict(sorted(symbol_bit_histogram.items())),
        "pointer_bit_histogram": dict(sorted(pointer_bit_histogram.items())),
        "total_bit_histogram": dict(sorted(total_bit_histogram.items())),
        "changed_symbol_cell_histogram": dict(sorted(changed_symbol_cell_histogram.items())),
        "average_symbol_bit_toggles": average_symbol_toggles,
        "average_pointer_bit_toggles": average_pointer_toggles,
        "average_total_bit_toggles": average_total_toggles,
        "max_symbol_bit_toggles": max(symbol_bit_histogram),
        "max_pointer_bit_toggles": max(pointer_bit_histogram),
        "max_total_bit_toggles": max(total_bit_histogram),
        "same_value_symbol_writes": changed_symbol_cell_histogram[0],
    }


def format_histogram(histogram: dict[int, int], num_states: int) -> str:
    return ", ".join(
        f"{toggles}:{count} ({100.0 * count / num_states:.1f}%)"
        for toggles, count in histogram.items()
    )


def print_summary(title: str, stats: dict[str, object]) -> None:
    print(title)
    for key, value in stats.items():
        if key.endswith("histogram"):
            continue
        print(f"  {key}: {value}")

    num_states = int(stats["num_states"])
    for key, value in stats.items():
        if key.endswith("histogram"):
            print(f"  {key}: {format_histogram(value, num_states)}")


def main() -> None:
    order = 4
    max_symbol = 8
    num_states = (max_symbol + 1) ** order

    binary_stats = binary_counter_switching_stats(num_states)
    onion_binary_pointer = onion_pointer_switching_stats(order, max_symbol, gray_pointer=False)
    onion_gray_pointer = onion_pointer_switching_stats(order, max_symbol, gray_pointer=True)

    print_summary(f"Binary modulo-{num_states} counter", binary_stats)
    print()
    print_summary(f"Onion counter with moving pointer (order={order}, max_symbol={max_symbol}, binary pointer)", onion_binary_pointer)
    print()
    print_summary(f"Onion counter with moving pointer (order={order}, max_symbol={max_symbol}, Gray pointer)", onion_gray_pointer)


if __name__ == "__main__":
    main()
