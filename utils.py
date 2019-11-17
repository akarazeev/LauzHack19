from Bio import SeqIO
from Bio import Align
import os
import numpy as np
from collections import OrderedDict


def load_patients(root_data='data/dataset_cp_pt_brain_dna/DNA', return_dict=False):
    patients = OrderedDict()
    for f in [f for f in os.listdir(root_data) if f.endswith('fastq')]:
        # print(f, end="\t... ")
        patients[f.split('.')[0]] = list(SeqIO.parse(os.path.join(root_data, f), "fastq"))
        # print("Done.")
        
    if return_dict:
        return patients
    else:
        patients_entries = np.array([[str(entry.seq)] + [k] for k in patients.keys() for entry in patients[k]])
        patients_identity = patients_entries[:, 1]
        patients_entries = patients_entries[:, 0]

        return patients_identity, patients_entries
    
    
def load_viruses(root_data='data/viral'):
    viruses = dict()
    for f in [f for f in os.listdir(root_data) if f.endswith('fna')]:
        # print(f, end="\t...")
        with open(os.path.join(root_data, f), 'r') as fileptr:
            text = fileptr.read().split('>')[1:]
            viruses.update({vir.split('\n')[0]: ''.join(vir.split('\n')[1:]) for vir in text})
            # print(" Done.")
    return viruses


def compute_virus_alignment_score(patients_entries, virus_dna, aligner):
    alignments_scores = np.zeros(patients_entries.shape[0])

    if len(virus_dna) != 0:
        for i, patient_dna in enumerate(patients_entries):
            if len(patient_dna) != 0:
                alignments = aligner.align(patient_dna, virus_dna)
                alignments_scores[i] = alignments.score
            else:
                alignments_scores[i] = 0

    return alignments_scores
