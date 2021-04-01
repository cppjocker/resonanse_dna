def calc_hydro_letter(two_letters):
    assert(len(two_letters) == 2)

    if two_letters[0] == 'N' or two_letters[1] == 'N':
        left_chain_start_code = 'N'
        return left_chain_start_code

    if two_letters == 'AA':
        left_chain_start_code = 'k'
    elif two_letters == 'AC':
        left_chain_start_code = 'm'
    elif two_letters == 'AG':
        left_chain_start_code = 'm'
    elif two_letters == 'AT':
        left_chain_start_code = 'f'
    elif two_letters == 'CA':
        left_chain_start_code = 'k'
    elif two_letters == 'CC':
        left_chain_start_code = 'm'
    elif two_letters == 'CG':
        left_chain_start_code = 'p'
    elif two_letters == 'CT':
        left_chain_start_code = 'm'
    elif two_letters == 'GA':
        left_chain_start_code = 'k'
    elif two_letters == 'GC':
        left_chain_start_code = 'm'
    elif two_letters == 'GG':
        left_chain_start_code = 'm'
    elif two_letters == 'GT':
        left_chain_start_code = 'm'
    elif two_letters == 'TA':
        left_chain_start_code = 'p'
    elif two_letters == 'TC':
        left_chain_start_code = 'k'
    elif two_letters == 'TG':
        left_chain_start_code = 'k'
    elif two_letters == 'TT':
        left_chain_start_code = 'k'
    else:
        print(two_letters)
        assert(False)

    return left_chain_start_code

def calc_hydro(seq, pos, nucl):
    left_start = seq[pos - 1] + nucl
    len_left = 0

    left_chain_start_letter = calc_hydro_letter(left_start)
    next_hydro_letter = left_chain_start_letter

    if next_hydro_letter in {'k', 'm'}:
        len_left = 1

    seq_orig = nucl
    seq_orig = seq[pos - 1].lower() + seq_orig
    seq_returned = next_hydro_letter.upper()

    for i in range(pos - 1, 0, -1):

        next_two_letters = seq[i - 1] + seq[i]
        next_hydro_letter = calc_hydro_letter(next_two_letters)

        seq_returned = next_hydro_letter + seq_returned
        seq_orig = seq[i - 1].lower() + seq_orig

        if next_hydro_letter not in {'k', 'm'}:
            break

        len_left = len_left + 1

    right_start = nucl + seq[pos + 1]
    len_right = 0

    right_chain_start_letter = calc_hydro_letter(right_start)
    next_hydro_letter = right_chain_start_letter

    if next_hydro_letter in {'k', 'm'}:
        len_right = 1


    seq_returned = seq_returned + next_hydro_letter.upper()
    seq_orig =  seq_orig + seq[pos + 1].lower()

    for i in range(pos + 1, len(seq)):
        next_two_letters = seq[i] + seq[i + 1]
        next_hydro_letter = calc_hydro_letter(next_two_letters)

        seq_returned = seq_returned + next_hydro_letter
        seq_orig =  seq_orig + seq[i + 1].lower()

        if next_hydro_letter not in {'k', 'm'}:
            break

        len_right = len_right + 1

    if left_chain_start_letter not in {'k', 'm'} or right_chain_start_letter not in {'k', 'm'}:
        total_len = max(len_left, len_right)
    else:
        total_len = len_left + len_right

    return total_len, seq_orig, seq_returned