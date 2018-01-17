import sys
import os
import re


def scan_motifs(motifs, genes, genes_compliment, max_length):
    """ Scan for promoters before genes.
    walks troh the motifs ands trys to find genes that are between de the max en minimum bp after the found promoter.
    When a gene is found it will me appended to a list with the promoter that fits.
    After all the motifs are looped the list with info will be returned.

    :param motifs: TYPE(list);  list with motif info from FIMO.
    :param genes: TYPE(List); List with info from the genes out the gbff file. On the plus streng.
    :param genes_compliment: TYPE(List); List with info from the genes out the gbff file. On the minus streng.
    :param max_length: TYPE(int); The max length before genes for the promoter.
    :return:
    """
    found_promoters = []
    for motif in motifs:
        if motif[5] == '+':
            for i in genes:
                minloc = int(re.sub('[^0-9]', '', i[0]
                                    .split('/')[0].split('..')[0])) - max_length
                maxloc = int(re.sub('[^0-9]', '', i[0]
                                    .split('/')[0].split('..')[0]))
                if int(motif[3]) >= minloc and int(motif[4]) < maxloc:
                    found_promoters.append((i, motif))
        if motif[5] == '-':
            for i in genes_compliment:
                minloc = int(re.sub('[^0-9]', '', i[0].split('/')[0].
                                    lstrip('complement()<>').split('..')[1]))
                maxloc = int(
                    re.sub('[^0-9]', '', i[0].split('/')[0]
                           .lstrip('complement()<>').split('..')[1])) + max_length
                if minloc < int(motif[3]) and int(motif[4]) <= maxloc:
                    found_promoters.append((i, motif))
    return found_promoters


def gbfiles(input_gbfile):
    """This funcion opens the genbank(flat) file.
    It only extracts all the info from the gene that are needed.
    The info from the genes are split in two lists. one list of
    gene on the + streng and the - streng.

    :param input_gbfile: type(str) the genbank file with info
    :return: gene_norm type(list), gene_compliment type(list)
    """

    with open(input_gbfile, "r") as file:
        gbfile = file.readlines()
        for line in range(0, len(gbfile)):
            if 'ORIGIN' in gbfile[line]:
                end_features = line
            if 'FEATURES' in gbfile[line]:
                start_featurs = line
    gene_info = ''.join([gbfile[linesnum] for linesnum
                         in range(start_featurs, end_features)]) \
        .split('     gene            ')
    gene_info.pop(0)
    gene_compliment, gene_norm = split_strengs(gene_info)
    return gene_norm, gene_compliment


def split_strengs(gene_info):
    """Here are the gene split in two groups. The + and - streng.
     If there are other like < > these gene are not counted
      cause the location is not know.


    :param gene_info: type(list); list with all the genes
    :return: gene_compliment type(list), gene_norm type(list)
    """
    gene_compliment = []
    gene_norm = []

    for i in gene_info:
        if ''.join(i.split())[0].isdigit():
            gene_norm.append([''.join(i.split())])
        elif ''.join(i.split())[0] == '<' or ''.join(i.split())[0] == '>':
            pass
        elif ''.join(i.split())[0] == 'c':
            gene_compliment.append([''.join(i.split())])
    return gene_compliment, gene_norm


def fil_file(promoters, distance, setting):
    """This function will give outpit to a file. In this file all the
     found gene after a promoter will be logged.
    There are two ways the output is written. easy and small. or extensive.
    Above each output there is a small explanation witch value which is.


    :param promoters: type(list)
    :param distance: type(int)
    :param setting: type(int)
    :return: None
    """
    with open('found_genes.txt', 'w+')as file:
        file.write(
            'Fimo promoter tool V 0.1\nTotal found genes: {}\nSetting:'
            ' {} bp before gene. output type: {}\n'.format(
                len(promoters), distance, setting))
        if setting == 0:
            file.write(
                'info:\n # motif_id	motif_alt_id	sequence_name	start	stop'
                '	strand	score	p-value'
                '	q-value	matched_sequence\n\n')
            for i in promoters:
                i[1].pop(1)
                file.write('Gene:\n' + ''.join(['\t' + i +
                                                '\n' for i in i[0][0]
                                               .replace('CDS', '/CDS\t')
                                               .split('/')]))
                file.write('\nMotif:\n' + ''.join(['\t' + i for i in i[1]]))
                file.write('+-----------------------------------------------+\n')
        else:
            file.write('\nOutput:\tgene_loc\tproduct\tid\tstart_motif'
                       '\tend_motif\tstreng\tmotif\n')
            for i in promoters:
                gen = [i for i in i[0][0].replace('CDS', '/CDS').split('/')]
                file.write(gen[0] + '\t' + gen[len(gen) - 2] + '\t'
                           + gen[len(gen) - 1] + '\t'
                           + i[1][3] + '\t' + i[1][
                               4] + '\t' + i[1][5] + '\t' + i[1][9])


def bad_files(error):
    """ This function will stop the program when the files are not present.
    The message print ands says witch file is not found.


    :param error: type(str); The name of the error of file.
    :return: None
    """
    print('One or more files are missing.')
    print('Blame {}'.format(error))
    exit(1)


def get_input():
    """This function will get al inputs form the terminal
     and will check if they are valid. if some error occurs
      in files bad_file() will be called and programs stops.
    When bp_before_gene or out_set are not a valid type there
     is a default value that will be used.

    :return: input_fimo type(str), input_gbfile type(str),
             bp_before_gene type(int), out_set type(int)
    """
    input_fimo = sys.argv[1]
    input_gbfile = sys.argv[2]
    try:
        bp_before_gene = int(sys.argv[3])
    except IndexError or TypeError:
        bp_before_gene = 300
    if not os.path.isfile(input_gbfile):
        bad_files(input_gbfile)
    if not os.path.isfile(input_fimo):
        bad_files(input_fimo)
    try:
        out_set = sys.argv[4]
    except IndexError or TypeError:
        out_set = 0

    return input_fimo, input_gbfile, bp_before_gene, out_set


def main():
    """This is the main function that runs everything.
    First the input is checked is its useful. When that's done there is
     feedback to the user what is happening. The FIMO file is opened and the
    input is written to a list. Now the same is done for the genes.
     When that's done the genes after a motif are found.
    Than the output is written to the file.

    :return: None
    """
    input_fimo, input_gbfile, bp_before_gene, settings = get_input()

    motifs = []

    print('Find FIMO motif hits V0.1\nbp before gene: {}\nOutput: {}\n\n'.format(bp_before_gene, settings))
    print('Opening motifs')
    with open(input_fimo, 'r') as fimo:
        for i in fimo.readlines():
            if i[0] != '#':
                motifs.append(i.split('\t'))

    print('[====>\t\t\t] 30 %\n\n\nOpening genbank file, extracting genes')

    gene, gene_compliment = gbfiles(input_gbfile)

    print('Finding genes\n[========>\t\t] 60 %\n\n')

    promoters = scan_motifs(motifs, gene, gene_compliment, bp_before_gene)

    fil_file(promoters, bp_before_gene, settings)

    print('Writing output to file\n[===============] 100% done\nOutput in found_gene.txt')


main()
