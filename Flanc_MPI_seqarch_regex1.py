from mpi4py import MPI
import argparse as ag
from scipy.spatial.distance import hamming
from Bio import SeqIO
import sys
import re

seq_start = 'gtcggggtttgtaccgtacaccact'
seq_end = 'ggtggttgaccagacaaaccacgac'
key = 'attP_Bxb1_reverse_complement_sequence'


def ag_pars():
    parser = ag.ArgumentParser()

    parser.add_argument('-g',
                        '--genome-input-file',
                        required=True,
                        dest='path_to_genome',
                        type=str,
                        help='Input file with gens')

    parser.add_argument('-o',
                        '--output-file',
                        required=True,
                        dest='new_file_name',
                        type=str,
                        help='Relative path to new folder')

    args = vars(parser.parse_args())
    return args


def fasta_read(attps_path: str) -> dict:
    # считывает файл fasta
    all_genes = dict()
    fasta_sequences = SeqIO.parse(attps_path, 'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        if name and sequence:
            all_genes[name] = sequence.upper()
        else:
            continue

    return all_genes


def repeat_search(genome: str, genome_header: str, regex: str, regex_re: str, key:str):
    # Сравнивает все кмеры генома с задаными attps и возвращает список
    data = []

    # [[genome_header, attp_header, str(-1)], seq, attp]
    seqs_ = re.findall(regex, genome, flags=0)
    seqs_r = re.findall(regex_re, genome, flags=0)
    
    for i in seqs_:
        data.append([[genome_header, key, str(1)], i, regex])
    
    for i in seqs_r:
        data.append([[genome_header, key, str(-1)], i, regex])
        
    return data


def repeat_process_mpi(gen_headers: list, gens: dict, new_file_name: str,
                       regex: str, regex_re: str, key:str):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    kmer_indx = rank
    gen_headers = tuple(gen_headers)

    with open(str(rank) + new_file_name, 'a+') as new_file:
        while kmer_indx < len(gen_headers):
            try:
                header = gen_headers[kmer_indx]
                genome = gens[header]
                print(kmer_indx)
                kmers = repeat_search(genome_header=header,
                                      genome=genome,
                                      regex=regex,
                                      regex_re=regex_re,
                                      key=key)
                if kmers:
                    pass
                else:
                    continue

                for i in kmers:
                    if i:
                        new_file.write(' '.join(['>'] + i[0] + ['\n']))
                        new_file.write(' '.join([i[1], '/', i[2], '\n']))
                    else:
                        pass
            except:
                print('Wrong way')
                pass
            finally:
                kmer_indx += size


def main():

    args = ag_pars()

    seq_re = seq_start + r'\w{2,15}' + seq_end
    seq_re_r = seq_end + r'\w{2,15}' + seq_start

    gens = fasta_read(args['path_to_genome'])
    gen_headers = tuple(gens.keys())

    repeat_process_mpi(gen_headers=gen_headers,
                       gens=gens,
                       new_file_name=args['new_file_name'],
                       regex=seq_re,
                       regex_re=seq_re_r,
                       key=key)


if __name__ == '__main__':
    main()