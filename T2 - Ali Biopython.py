from Bio import Align
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices

def sw_biopython(self, seq1, seq2):
    alineador_f = PairwiseAligner() # funcion alineador comprada
    alineador_f.mode = "local"
    alineador_f.substitution_matrix = substitution_matrices.load("Matrices_name")

    alineador_f.open_gap_score = self.config["gap_penalty"] # en el config
    alineador_f.extend_gap_score = self.config["gap_penalty"]# en el config

    # accion PairwiseAligner().align() PERO TUNEADO POR ESO NO DIREECTO
    alineamientos_acc = alineador_f.align(seq1,seq2)
    return alineamientos_acc[0].score
