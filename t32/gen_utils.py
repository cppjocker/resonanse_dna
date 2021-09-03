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
