from Bio import SeqIO
from Bio import Align
import os


def load_patients(root_data='data/dataset_cp_pt_brain_dna/DNA'):
    patients = dict()
    for f in [f for f in os.listdir(root_data) if f.endswith('fastq')]:
        print(f, end="\t... ")
        patients[f.split('.')[0]] = list(SeqIO.parse(os.path.join(root_data, f), "fastq"))
        print("Done.")
    return patients
    
    
def load_viruses(root_data='data/viral'):
    viruses = dict()
    for f in [f for f in os.listdir(root_data) if f.endswith('fna')]:
        print(f, end="\t...")
        with open(os.path.join(root_data, f), 'r') as fileptr:
            text = fileptr.read().split('>')[1:]
            viruses.update({vir.split('\n')[0]: ''.join(vir.split('\n')[1:]) for vir in text})
            print(" Done.")
    return viruses

