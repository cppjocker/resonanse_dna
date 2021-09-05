


def calc_AT_metric1(seq, pos, nucl):
    if not ( nucl == 'A' and seq[pos + 1] == 'T' ):
        return -1, -1, -1

    l1 = 0
    for walk_left in range(pos - 2, 0, -1):
        if (seq[walk_left] == 'A' and seq[walk_left + 1] == 'T'):
            l1 = (pos - 1) - (walk_left + 1)
            break

    l2 = 0
    for walk_right in range (pos + 1, len(seq) ):
        if ( seq[walk_right] == 'A' and seq[walk_right + 1] == 'T' ):
            l2 = (walk_right - 1) - (pos + 1)
            break


    return l1, l2, l1*l2

def calc_AT_metric2(seq, pos, nucl):
    if not ( seq[pos - 1] == 'A' and nucl  == 'T' ):
        return -1, -1, -1

    l1 = 0
    for walk_left in range(pos - 3, 0, -1):
        if (seq[walk_left] == 'A' and seq[walk_left + 1] == 'T'):
            l1 = (pos - 2) - (walk_left + 1)
            break

    l2 = 0
    for walk_right in range (pos + 1, len(seq) ):
        if ( seq[walk_right] == 'A' and seq[walk_right + 1] == 'T' ):
            l2 = (walk_right - 1) - (pos)
            break


    return l1, l2, l1*l2

def calc_AT_metric(seq, pos, nucl, code):
    if code == 1:
        return calc_AT_metric1(seq, pos, nucl)
    elif code == 2:
        return calc_AT_metric2(seq, pos, nucl)


def calc_shoulder_metric(seq, pos, nucl, code):
    if not ( nucl == code[0] and seq[pos + 1] == code[1] ):
        return -1, -1, -1, False

    l1 = 0
    for walk_left in range(pos - 2, 0, -1):
        if (seq[walk_left] == code[0] and seq[walk_left + 1] == code[1]):
            l1 = (pos - 1) - (walk_left + 1)
            break

    l2 = 0
    for walk_right in range (pos + 1, len(seq) ):
        if ( seq[walk_right] == code[0] and seq[walk_right + 1] == code[1] ):
            l2 = (walk_right - 1) - (pos + 1)
            break

    AT_in = False

    if "AT" in seq[walk_left:(walk_right + 1)]:
        AT_in = True


    return l1, l2, l1*l2, AT_in