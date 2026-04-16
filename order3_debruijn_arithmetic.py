Order3State = tuple[int, int, int]

def integer_cube_root(x: int) -> int:
    m = int(x ** (1 / 3))
    if m ** 3 <= x < (m + 1) ** 3:
        return m
    else:
        return m + 1


def parse_state(raw_state: Order3State | str | int) -> Order3State:
    if isinstance(raw_state, int):
        raw_state = f'{raw_state:03}'
    if isinstance(raw_state, str):
        state = int(raw_state[0], 36), int(raw_state[1], 36), int(raw_state[2], 36)
    else:
        state = raw_state
    return state


def rho3(state: Order3State | str | int) -> int:
    a, b, c = parse_state(state)
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
    raise ValueError(f'Unable to calculate index of {state}')


def inverse_rho3(index: int) -> Order3State:
    layer = integer_cube_root(index)
    full_offset = index - layer ** 3
    word_index, offset = divmod(full_offset, 3)
    primary_word_index, secondary_word_index = divmod(word_index, layer) if layer else (0, 0)
    next_layer_index = (layer + 1) ** 3
    if index == next_layer_index - 3:
        return layer, layer, layer
    elif index == next_layer_index - 2:
        return layer, layer, 0
    elif index == next_layer_index - 1:
        return layer, 0, 0
    elif offset == 0:
        return secondary_word_index, primary_word_index, layer
    elif offset == 1:
        return primary_word_index, layer, (secondary_word_index + 1) % layer
    elif offset == 2:
        if (secondary_word_index + 1) % layer == 0:
            return layer, 0, primary_word_index + 1
        else:
            return layer, (secondary_word_index + 1) % layer, primary_word_index
    raise ValueError(f'Unable to calculate state-3 at {index}')


def lambda3(state: Order3State | str | int) -> tuple[int, int, int]:
    a, b, c = parse_state(state)
    if  b == c == 0:
        return a + 1, 0, -1
    elif a == b and c == 0:
        return a, a ** 2 + a, -1
    elif a == b == c:
        return a, a ** 2 + a - 1, 1
    elif c >= b and c > a:  # c is the pivot
        return c, b * c + a, 0
    elif b >= a and b > c:  # b is the pivot
        return b, a * b + (c - 1) % b, 1
    elif a >= c and a > b:  # a is the pivot
        return a, c * a + b, -1
    raise ValueError(f'Cannot calculate lambda3({state})')


def inverse_lambda3(layer: int, word_index: int, offset: int) -> Order3State:
    if offset == 0:
        return word_index % layer, word_index // layer, layer
    elif offset == 1:
        if word_index == (layer ** 2 + layer - 1):
            return layer, layer, layer
        else:
            return word_index // layer, layer, (word_index + 1) % layer
    else:  # offset == -1
        if word_index == 0:
            return layer - 1, 0, 0
        elif word_index == (layer ** 2 + layer):
            return layer, layer, 0
        else:
            return layer, word_index % layer, word_index // layer


def carry_down(layer: int, word_index: int, offset: int) -> tuple[int, int, int]:
    while word_index < 0:
        layer, word_index, offset = layer - 1, word_index + layer ** 2 - layer + (offset + 2) // 3, (offset - 1) % 3 - 1
    return layer, word_index, offset


def add_states3(state1: Order3State | int | str, state2: Order3State | int | str) -> Order3State:
    layer1, word_index1, offset1 = lambda3(state1)
    layer2, word_index2, offset2 = lambda3(state2)
    combined_layer = layer1 + layer2
    combined_word_index = word_index1 + word_index2 - layer1 * layer2 * combined_layer + (offset1 + offset2 + 1) // 3
    combined_offset = (offset1 + offset2 + 1) % 3 - 1
    combined_layer, combined_word_index, combined_offset = carry_down(combined_layer, combined_word_index, combined_offset)
    return inverse_lambda3(combined_layer, combined_word_index, combined_offset)


def carry_up(layer: int, word_index: int, offset: int) -> tuple[int, int, int]:
    while word_index > layer ** 2 + layer - offset // 2 - 1:
        layer, word_index, offset = layer + 1, word_index - layer ** 2 - layer - (2 - offset) // 3, offset % 3 - 1
    return layer, word_index, offset


def mul_states3(state1: Order3State | int | str, state2: Order3State | int | str) -> Order3State:
    layer1, word_index1, offset1 = lambda3(state1)
    layer2, word_index2, offset2 = lambda3(state2)
    combined_layer = layer1 * layer2
    offset_remainder = offset2 * layer1 ** 3 + offset1 * layer2 ** 3 + offset1 * offset2
    combined_word_index = layer1 ** 3 * word_index2 + layer2 ** 3 * word_index1 + 3 * word_index1 * word_index2 + word_index1 * offset2 + word_index2 * offset1 + (offset_remainder - 2) // 3 + 1
    combined_offset = (offset_remainder + 1) % 3 - 1
    if combined_word_index > 0:
        combined_layer, combined_word_index, combined_offset = carry_up(combined_layer, combined_word_index, combined_offset)
    else:
        combined_layer, combined_word_index, combined_offset = carry_down(combined_layer, combined_word_index, combined_offset)
    return inverse_lambda3(combined_layer, combined_word_index, combined_offset)


def divmod_states3(state1: Order3State, state2: Order3State) -> tuple[Order3State, Order3State]:
    divisor = rho3(state2)
    if divisor == 0:
        raise ZeroDivisionError("division by zero onion state")
    quotient, remainder = divmod(rho3(state1), divisor)
    return inverse_rho3(quotient), inverse_rho3(remainder)


def floordiv_states3(state1: Order3State, state2: Order3State) -> Order3State:
    return divmod_states3(state1, state2)[0]


def mod_states3(state1: Order3State, state2: Order3State) -> Order3State:
    return divmod_states3(state1, state2)[1]
