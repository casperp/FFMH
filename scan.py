import sys
import re

"""Find fimo motif hits
FFFH
"""


def scan_motifs(motifs, genes, genes_compliment, max_lenght):
    found_promoters = []
    for motif in motifs:
        if motif[5] == '+':
            for i in genes:
                minloc = int(re.sub('[^0-9]', '',i[0].split('/')[0].split('..')[0])) - max_lenght
                maxloc = int(re.sub('[^0-9]', '',i[0].split('/')[0].split('..')[0]))
                if int(motif[3]) >= minloc and int(motif[4]) < maxloc:
                    found_promoters.append((i, motif))
        if motif[5] == '-':
            for i in genes_compliment:
                minloc = int(re.sub('[^0-9]', '', i[0].split('/')[0].lstrip('complement()<>').split('..')[1]))
                maxloc = int(
                    re.sub('[^0-9]', '', i[0].split('/')[0].lstrip('complement()<>').split('..')[1])) + max_lenght
                if minloc < int(motif[3]) and int(motif[4]) <= maxloc:
                    found_promoters.append((i, motif))
    return found_promoters


def gbfiles(input_gbfile):
    gene_compliment = []
    gene_norm = []
    with open(input_gbfile, "r") as file:
        gbfile = file.readlines()
        for line in range(0, len(gbfile)):
            if 'ORIGIN' in gbfile[line]:
                end_features = line
            if 'FEATURES' in gbfile[line]:
                start_featurs = line
    gene_info = ''.join([gbfile[linesnum] for linesnum
                         in range(start_featurs, end_features)]).split('     gene            ')
    gene_info.pop(0)
    for i in gene_info:
        if ''.join(i.split())[0].isdigit():
            gene_norm.append([''.join(i.split())])
        elif ''.join(i.split())[0] == '<' or ''.join(i.split())[0] == '>':
            pass
        elif ''.join(i.split())[0] == 'c':
            gene_compliment.append([''.join(i.split())])
    return gene_norm, gene_compliment


def fil_file(promoters,settings):
    with open('found_genes.txt', 'w+')as file:
        file.write('Fimo promoter tool V 0.1\nTotal found genes: {}\n\n'.format(len(promoters)))
        file.write(
            'info:\n # motif_id	motif_alt_id	sequence_name	start	stop	strand	score	p-value	q-value	matched_sequence\n\n')
        for i in promoters:
            i[1].pop(1)
            # print([i for i in i[0][0].replace('CDS', '/CDS\t').split('/')])
            file.write('Gene:\n' + ''.join(['\t' + i + '\n' for i in i[0][0].replace('CDS', '/CDS\t').split('/')]))
            file.write('\nMotif:\n' + ''.join(['\t' + i for i in i[1]]))
            file.write('+-----------------------------------------------+\n')

    pass


def get_input():
    input_fimo = sys.argv[1]
    input_gbfile = sys.argv[2]
    try:
        bp_before_gene = int(sys.argv[3])
    except TypeError:
        bp_before_gene = 300
    return input_fimo, input_gbfile, bp_before_gene


def main():
    input_fimo,input_gbfile, bp_before_gene = get_input()
    motifs = []
    with open(input_fimo, 'r') as fimo:
        for i in fimo.readlines():
            if i[0] != '#':
                motifs.append(i.split('\t'))
    gene, gene_compliment = gbfiles(input_gbfile)

    promoters = scan_motifs(motifs, gene, gene_compliment, bp_before_gene)
    fil_file(promoters, bp_before_gene)


main()
