


def calc_AT_metric(seq, pos, nucl):
    if not ( nucl == 'A' and seq[pos + 1] == 'T' ):
        return -1

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


    return l1 * l2