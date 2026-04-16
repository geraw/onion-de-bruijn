Order2State = tuple[int, int]

def parse_state(raw_state: Order2State | str | int) -> Order2State:
    if isinstance(raw_state, int):
        raw_state = f'{raw_state:02}'
    if isinstance(raw_state, str):
        state = int(raw_state[0], 36), int(raw_state[1], 36)
    else:
        state = raw_state
    return state


def rho2(state: Order2State | str | int) -> int:
    a, b = parse_state(state)
    if a < b:
        return b ** 2 + 2 * a
    elif 0 < b <= a:
        return a ** 2 + 2 * b - 1
    else:
        return a ** 2 + 2 * a


def inverse_rho2(index: int) -> Order2State:
    layer = int(index ** 0.5)
    word_index = index - layer ** 2
    if word_index % 2:
        return layer, (word_index + 1) // 2
    elif word_index < 2 * layer:
        return word_index // 2, layer
    else:
        return layer, 0


def lambda2(state: Order2State | str | int) -> tuple[int, int, int]:
    a, b = parse_state(state)
    if a < b:
        return b, a, 0
    elif b == 0:
        return a, a, 0
    else:
        return a, b, 1


def inverse_lambda2(layer: int, word_index: int, offset: int) -> Order2State:
    if offset == 0:
        if word_index < layer:
            return word_index, layer
        else:
            return layer, 0
    else:
        return layer, word_index


def carry_down(layer: int, word_index: int, offset: int) -> tuple[int, int, int]:
    while (word_index + 1 - offset) <= 0:
        layer, word_index, offset = layer - 1, word_index + layer - offset, 1 - offset
    return layer, word_index, offset


def carry_up(layer: int, word_index: int, offset: int) -> tuple[int, int, int]:
    while word_index > layer:
        layer, word_index, offset = layer + 1, word_index - layer - offset, 1 - offset
    return layer, word_index, offset


def add_states2(state1: Order2State | str | int, state2: Order2State | str | int) -> Order2State:
    layer1, word_index1, offset1 = lambda2(state1)
    layer2, word_index2, offset2 = lambda2(state2)
    combined_layer = layer1 + layer2
    combined_word_index = word_index1 + word_index2 - layer1 * layer2 - (offset1 + offset2) // 2
    combined_offset = (offset1 + offset2) % 2
    combined_layer, combined_word_index, combined_offset = carry_down(combined_layer, combined_word_index, combined_offset)
    return inverse_lambda2(combined_layer, combined_word_index, combined_offset)


def mul_states2(state1: Order2State | str | int, state2: Order2State | str | int) -> Order2State:
    layer1, word_index1, offset1 = lambda2(state1)
    layer2, word_index2, offset2 = lambda2(state2)
    combined_layer = layer1 * layer2
    offset_remainder = offset2 * layer1 ** 2 + offset1 * layer2 ** 2 - offset1 * offset2
    combined_word_index = layer1 ** 2 * word_index2 + layer2 ** 2 * word_index1 + 2 * word_index1 * word_index2 - word_index1 * offset2 - word_index2 * offset1 - offset_remainder // 2
    combined_offset = offset_remainder % 2
    combined_layer, combined_word_index, combined_offset = carry_up(combined_layer, combined_word_index, combined_offset)
    return inverse_lambda2(combined_layer, combined_word_index, combined_offset)


def divmod_states2(state1: Order2State, state2: Order2State) -> tuple[Order2State, Order2State]:
    divisor = rho2(state2)
    if divisor == 0:
        raise ZeroDivisionError("division by zero onion state")
    quotient, remainder = divmod(rho2(state1), divisor)
    return inverse_rho2(quotient), inverse_rho2(remainder)


def floordiv_states2(state1: Order2State, state2: Order2State) -> Order2State:
    return divmod_states2(state1, state2)[0]


def mod_states2(state1: Order2State, state2: Order2State) -> Order2State:
    return divmod_states2(state1, state2)[1]
