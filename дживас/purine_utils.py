def calc_RY(seq, pos, nucl):
    R_chain_len = 0
    Y_chain_len = 0

    ry_nucl = RY(nucl)


    if ry_nucl == 'R':
        R_chain_len = 1
        Y_chain_len = 0
        ry_nucl_compl = 'Y'
    else:
        R_chain_len = 0
        Y_chain_len = 1

        ry_nucl_compl = 'R'

    chain_len = 0

    for i in range(pos - 1, 0, -1):
        if RY(seq[i]) != ry_nucl:
            break

        chain_len+=1

    for i in range(pos + 1, len(seq)):
        if RY(seq[i]) != ry_nucl:
            break

        chain_len+=1

    if ry_nucl == 'R':
        R_chain_len += chain_len
    else:
        Y_chain_len += chain_len

    chain_len_compl_neg = 0
    chain_len_compl_pos = 0

    for i in range(pos - 1, 0, -1):
        if RY(seq[i]) != ry_nucl_compl:
            break

        chain_len_compl_neg+=1

    for i in range(pos + 1, len(seq)):
        if RY(seq[i]) != ry_nucl_compl:
            break

        chain_len_compl_pos+=1

    chain_len_compl = max(chain_len_compl_neg, chain_len_compl_pos)

    if ry_nucl_compl == 'R':
        R_chain_len += chain_len_compl
    else:
        Y_chain_len += chain_len_compl

    return R_chain_len, Y_chain_len

def RY(nucl):
    if (nucl == 'A') or (nucl == 'G'):
        return 'R'
    else:
        return 'Y'