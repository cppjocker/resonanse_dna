def calc_cons_scores(pos_SNP, f, last_pos, step):
    score = -1.0

    while last_pos < pos_SNP:

        line = f.readline()

        if not line:
            return False, last_pos, step, -1

        last_pos = last_pos + step

        parts = line.split()
        assert( len(parts) == 4 or len(parts) == 1 )

        if len(parts) == 4:
            assert( parts[0] == 'fixedStep')

            chrom = parts[1]
            start = parts[2]
            new_step = parts[3]

            start_spl = start.split('=')
            new_step_spl = new_step.split('=')

            assert(chrom.split('=')[0] == 'chrom')
            assert(start_spl[0] == 'start')
            assert(new_step_spl[0] == 'step')

            last_pos = int( start_spl[1] )
            step = int(new_step_spl[1])

            line = f.readline()
            assert (line)

            parts = line.split()
            assert (len(parts) == 1)

        assert(len(parts) == 1)

        score = float(line)

    if last_pos > pos_SNP:
        return False, last_pos, step, -1

    return True, last_pos, step, score





def calc_gen_markup(pos, idx_trans, df_gtf_cur, df_gtf_cur_trans):
    found_intersec = False
    return_pos = 0

    has_exon = False
    has_intron = 0
    has_3utr = 0
    has_5utr = 0

    has_itrgen = False

    idx_trans_start = idx_trans

    while True:
        if idx_trans >= df_gtf_cur_trans.shape[0] - 1:
            break

        start_trans = df_gtf_cur_trans.iat[idx_trans, 3]
        end_trans = df_gtf_cur_trans.iat[idx_trans, 4]

        next_trans_start = df_gtf_cur_trans.iat[idx_trans + 1, 3]

        assert (start_trans <= next_trans_start)

        if start_trans > pos:
            break

        if (not found_intersec) and (start_trans <= pos <= end_trans):
            found_intersec = True
            return_pos = idx_trans

        if start_trans <= pos <= end_trans:
            orig_index_start = df_gtf_cur_trans.index[idx_trans]
            orig_index_end = df_gtf_cur_trans.index[idx_trans + 1]

            has_exon_cur = False
            has_3utr_cur = 0
            has_5utr_cur = 0

            for i in range(orig_index_start, orig_index_end):
                cur_feature = df_gtf_cur.iat[i, 2]
                cur_feature_start = df_gtf_cur.iat[i, 3]
                cur_feature_end = df_gtf_cur.iat[i, 4]

                if not (cur_feature in {'CDS', '3UTR', '5UTR'}):
                    continue

                if not (cur_feature_start <= pos <= cur_feature_end):
                    continue

                if cur_feature == 'CDS':
                    has_exon_cur = True

                if cur_feature == '3UTR':
                    has_3utr_cur = 1

                if cur_feature == '5UTR':
                    has_5utr_cur = 1

            has_exon = has_exon or has_exon_cur

            if has_3utr_cur == 1:
                has_3utr = has_3utr_cur

            if has_5utr_cur == 1:
                has_5utr = has_5utr_cur

            if (not has_exon_cur) and (has_3utr_cur + has_5utr_cur == 0):
                has_intron = True

        idx_trans = idx_trans + 1

    idx_trans_end = idx_trans

    idx_trans_start = max(idx_trans_start - 20, 0)
    idx_trans_end = min(idx_trans_end + 20, df_gtf_cur_trans.shape[0])

    dist = 1000000000000
    strand = 'xx'

    for i in range(idx_trans_start, idx_trans_end):
        if df_gtf_cur_trans.iat[i, 6] == '+':  # strand
            cur_dist = abs(pos - df_gtf_cur_trans.iat[i, 3])  # start
            if pos < df_gtf_cur_trans.iat[i, 3]:
                cur_strand = 'up'
            else:
                cur_strand = 'down'
        elif df_gtf_cur_trans.iat[i, 6] == '-':
            cur_dist = abs(pos - df_gtf_cur_trans.iat[i, 4])  # end
            if pos > df_gtf_cur_trans.iat[i, 3]:
                cur_strand = 'up'
            else:
                cur_strand = 'down'
        else:
            assert (False)

        if cur_dist < dist:
            dist = cur_dist
            strand = cur_strand

    if found_intersec:
        idx_trans = return_pos

    weird_sum = has_intron + has_3utr + has_5utr

    if (not has_exon) and (weird_sum == 0):
        has_itrgen = True

    gen = 'unk'

    if has_itrgen:
        gen = 'igc'
    elif has_exon:
        gen = 'CDS'
    elif has_intron and weird_sum == 1:
        gen = 'tro'
    elif has_3utr and weird_sum == 1:
        gen = 'tpu'
    elif has_5utr and weird_sum == 1:
        gen = 'fpu'

    return idx_trans, gen, dist, strand
