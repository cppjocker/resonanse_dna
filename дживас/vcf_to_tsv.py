import os

src_dir_path = './input_cacao'
trg_dir_path = './output_cacao'

if not os.path.exists(trg_dir_path):
    os.makedirs(trg_dir_path)

src_dir_names = os.listdir(src_dir_path)
for dir in src_dir_names:
    src_file_names = os.listdir(src_dir_path + '/' + dir)
    trg_file_name = dir + '_tab.tsv'
    val_lines = 'CHROM' + '\t' + 'POS' + '\t' + 'N_ALLELES' + '\t' + 'A1' + '\t' + 'A1FREQ' + '\t' + 'A2' + '\t' + 'A2FREQ' + '\t' + 'REF' + '\n'
    for src_file_name in src_file_names:
        with open(os.path.join(src_dir_path,dir,src_file_name)) as src_file_opened:
            for line in src_file_opened:
                print(src_file_name, line)
                if line.startswith('#') == False:
                    chrom = line.split()[0]
                    pos = line.split()[1]
                    ref = line.split()[3]
                    a1 = line.split()[4]
                    n_al = 1
                    a2 = '.'
                    af1 = '.'
                    af2 = '.'
                    if ',' in line.split()[4]:
                        a1 = line.split()[4].split(',')[0]
                        a2 = line.split()[4].split(',')[1]
                        n_al = len(line.split()[4].split(','))
                    for i in line.split()[7].split(';'):
                        if 'AF' in i:
                            af = i.split('=')[1]
                            af1 = af.split(',')[0]
                            if len(af.split(','))>1:
                                af2 = af.split(',')[1]
                    val_lines += str(chrom) + '\t' + str(pos)+ '\t' + str(n_al) + '\t' + str(a1) + '\t' + str(af1) + '\t' + str(a2) + '\t' + str(af2) + '\t' + str(ref) + '\n'
    with open(os.path.join(trg_dir_path, trg_file_name), 'w') as trg_file_opened:
        trg_file_opened.write(val_lines)



