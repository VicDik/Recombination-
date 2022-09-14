from mpi4py import MPI
import argparse as ag
from scipy.spatial.distance import hamming
from Bio import SeqIO
import sys


def ag_pars():
    parser = ag.ArgumentParser()

    parser.add_argument(
        '-a',
        '--attp-input-file',
        required=True,
        dest='path_to_attps',
        type=str,
        help='Input file with attps')
    
    parser.add_argument(
        '-g',
        '--genome-input-file',
        required=True,
        dest='path_to_genome',
        type=str,
        help='Input file with genome')

    parser.add_argument(
        '-s',
        '--similarity',
        required=True,
        dest='similarity',
        type=float,
        help='procent of similarity for attps and k-mer')

    parser.add_argument(
        '-o',
        '--output-file',
        required=True,
        dest='new_file_name',
        type=str,
        help='Relative path to new folder')

    args = vars(parser.parse_args())
    return args

def fasta_read(attps_path:str) -> dict:
    # считывает файл fasta
    all_genes = dict()
    fasta_sequences = SeqIO.parse(attps_path,'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        if name and sequence:
            all_genes[name] = sequence.upper()
        else:
            continue

    return all_genes

def kmers_sim_calcs(genome:str, attps:dict, genome_header:str, similarity=0.7, min_len_attp = 15):
    # Сравнивает все кмеры генома с задаными attps и возвращает список
    data = []
    for attp_header in attps:
        attp = attps[attp_header]
        k = len(attp)
        for i in range(len(genome) - k + 1):
            seq = genome[i: i+k]
            hamm_distance = 1 - hamming(list(attp), list(seq))
            hamm_distance_r = 1 - hamming(list(attp), list(seq[::-1]))
            if k < min_len_attp:
                break
            #elif attp == seq:
            elif hamm_distance >= similarity:
                data.append([[genome_header, attp_header, str(1)], seq, attp])
                continue
            #elif attp[::-1] == seq:
            elif hamm_distance_r >= similarity:
                data.append([[genome_header, attp_header, str(-1)], seq, attp])
            else:
                pass
    return data

def kmer_process_mpi(attps:dict, gene_headers:list, genes:dict, new_file_name: str, similarity:float):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    kmer_indx = rank
    gene_headers = tuple(gene_headers)
    
    with open(str(rank) + new_file_name, 'a+')  as new_file:
        while kmer_indx < len(gene_headers):
            try:
                header = gen_headers[kmer_indx]
                genome = gens[header]
                print(kmer_indx)
                kmers = kmers_sim_calcs(genome_header=header, genome=genome, attps=attps, similarity=similarity)
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
    
    attps = fasta_read(args['path_to_attps'])
    gens = fasta_read(args['path_to_genome'])
    gen_headers = tuple(gens.keys())
    
    kmer_process_mpi(attps=attps,
                     gen_headers=gen_headers,
                     gens=gens,
                     new_file_name=args['new_file_name'],
                     similarity=float(args['similarity']))
    
    

if __name__ == '__main__':
    main()
