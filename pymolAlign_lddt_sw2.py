import numpy as np
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.Polypeptide import PPBuilder
from pymol import cmd
import subprocess
import glob
import sys

def pdb_extract_distances(pdbfile):
    
    #read with Bio.python
    if pdbfile[-4:]==".pdb": parser = PDBParser(QUIET=True)
    if pdbfile[-4:]==".cif": parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", pdbfile)
    
    # extract info
    resid = [atom.get_parent().get_full_id()[3][1] for atom in structure.get_atoms() if atom.get_name() == "CA"] # stores residue number for each alpha carbon
    ca50_atoms = [atom for atom in structure.get_atoms() if atom.get_name() == "CA" and atom.get_bfactor() > 50] # stores the entries of alpha carbon atoms with bfactor (plddt)>50
    resid50 = [atom.get_parent().get_full_id()[3][1] for atom in structure.get_atoms() if atom.get_name() == "CA" and atom.get_bfactor() > 50] # stores residue number for each alpha carbon with bfactor(plddt) > 50
    aaseq = str(PPBuilder().build_peptides(structure)[0].get_sequence()) # stores one-letter code of all residues
    length = len(resid)

    # Find neighbors within the distance cutoff, as query structure, +4, as reference structure
    distance_cutoff = 15.0
    total_residues = len(ca50_atoms)
    allref_distances={}
    allquery_distances={}
    for i, atom in zip(resid50,ca50_atoms):
        refdistances = []
        querydistances = []
        for j, other_atom in zip(resid50,ca50_atoms):
            if i != j:  # Exclude self
                dist = np.linalg.norm(atom.coord - other_atom.coord)
                if dist < distance_cutoff + 4:
                    refdistances.append(dist)
                if dist < distance_cutoff:
                    querydistances.append(dist)

        allref_distances[i]=sorted(refdistances)
        allquery_distances[i]=sorted(querydistances)        
    

    return allref_distances, allquery_distances, length, aaseq, resid, resid50

def calculate_lddtmatchmatrix(s1_len,s2_len,s1_resid,s2_resid,s1_resid50, s2_resid50,s1_qdist,s2_rdist):
    
    mismatch_penalty = -1
    match_matrix = [[mismatch_penalty] * (s1_len+1) for _ in range(s2_len+1)]

    # Build match matrix with calculate lDDT scores averages for each threshold
    distance_cutoff=15.0
    thresholds=[0.5, 1, 2, 4]
    threshold_weight = [2.25,1,0.5,0.25] # so that closer atoms weight more
    for i in s1_resid50:
        for j in s2_resid50:
            lddt_threshold_average=0
            for w,threshold in enumerate(thresholds):
                matches = sum(abs(r - m) < threshold for r, m in zip(s2_rdist[j], s1_qdist[i]))
                lenlargestlist=len(s1_qdist[i])
                lencutoffrefdistances=sum(np.array(s2_rdist[j])<distance_cutoff)
                if lencutoffrefdistances > len(s1_qdist[i]): lenlargestlist=lencutoffrefdistances
                if lenlargestlist!=0:
                    lddt_score = matches / lenlargestlist
                else:
                    lddt_score = 1
                lddt_threshold_average += lddt_score * threshold_weight[w]
            lddt_threshold_average /= len(thresholds)
            match_matrix[s2_resid.index(j)][s1_resid.index(i)] = lddt_threshold_average
    
    return match_matrix
    
def smith_waterman_affine_nonoverlapping(seq1, seq2,match_matrix,min_score):
    gap_open=-1
    gap_extend=-0.5

    original_seq1 = list(seq1)
    original_seq2 = list(seq2)
    used_mask1 = [False] * len(seq1)
    used_mask2 = [False] * len(seq2)
    alignments = []

    while True:
        n = len(seq1)
        m = len(seq2)

        # Create score matrices
        M = np.zeros((n + 1, m + 1))
        Ix = np.full((n + 1, m + 1), -np.inf)
        Iy = np.full((n + 1, m + 1), -np.inf)
        traceback = np.zeros((n + 1, m + 1), dtype=int)

        NONE, DIAG, UP, LEFT = 0, 1, 2, 3
        max_score = 0
        max_pos = None

        # Fill matrices
        for i in range(1, n + 1):
            if seq1[i - 1] is None:
                continue
            for j in range(1, m + 1):
                if seq2[j - 1] is None:
                    continue

                Ix[i][j] = max(Ix[i - 1][j] + gap_extend, M[i - 1][j] + gap_open + gap_extend)
                Iy[i][j] = max(Iy[i][j - 1] + gap_extend, M[i][j - 1] + gap_open + gap_extend)

                match_score = M[i - 1][j - 1] + match_matrix[j-1][i-1]
                M[i][j] = max(0, match_score, Ix[i][j], Iy[i][j])

                if M[i][j] == match_score:
                    traceback[i][j] = DIAG
                elif M[i][j] == Ix[i][j]:
                    traceback[i][j] = UP
                elif M[i][j] == Iy[i][j]:
                    traceback[i][j] = LEFT
                else:
                    traceback[i][j] = NONE

                if M[i][j] > max_score:
                    max_score = M[i][j]
                    max_pos = (i, j)

        if max_score < min_score or max_pos is None:
            break

        # Traceback
        i, j = max_pos
        aligned1, aligned2 = [], []
        aligned1_resid, aligned2_resid = [], []
        path1_indices, path2_indices = set(), set()

        #print('New alignment')
        while M[i][j] > 0:
            if seq1[i - 1] is None or seq2[j - 1] is None:
                break

            path1_indices.add(i - 1)
            path2_indices.add(j - 1)

            if traceback[i][j] == DIAG:
                aligned1.insert(0, seq1[i - 1])
                aligned1_resid.insert(0,i-1)
                aligned2.insert(0, seq2[j - 1])
                aligned2_resid.insert(0,j-1)
                i, j = i - 1, j - 1
            else:
                break
            #print(seq1[i],seq2[j],M[i][j],match_matrix[j][i])


        if not aligned1 or not aligned2:
            break  # invalid or empty alignment

        # Store alignment
        alignments.append((aligned1, aligned1_resid, aligned2, aligned2_resid, max_score))

        # Mask out used indices in original sequences
        for idx in path1_indices:
            used_mask1[idx] = True
        for idx in path2_indices:
            used_mask2[idx] = True

        # Apply masks
        seq1 = [None if used_mask1[i] else original_seq1[i] for i in range(len(seq1))]
        seq2 = [None if used_mask2[j] else original_seq2[j] for j in range(len(seq2))]

    return alignments

def super_lddt(pdbs=[],min_score=50):
    
    if len(pdbs)!=2:
        sys.exit("Please enable only 2 objects")
    
    pdb1=glob.glob(pdbs[0]+'.*')[0]
    pdb2=glob.glob(pdbs[1]+'.*')[0]

    s1_refdistances, s1_querydistances, s1_length, s1_aaseq, s1_resid, s1_resid50 = pdb_extract_distances(pdb1)
    s2_refdistances, s2_querydistances, s2_length, s2_aaseq, s2_resid, s2_resid50 = pdb_extract_distances(pdb2)

    match_matrix = calculate_lddtmatchmatrix(s1_length,s2_length,s1_resid,s2_resid,s1_resid50,s2_resid50,s1_querydistances,s2_refdistances)
    print(len(match_matrix),len(match_matrix[0]))
    print(len(s1_aaseq),len(s2_aaseq))
    results = smith_waterman_affine_nonoverlapping(s1_aaseq, s2_aaseq, match_matrix, min_score)

    supers=pdbs[0]+'_supers'
    print("Creating a new object of %s, %s, which will be loaded with different alignments as different states." % (pdbs[0],supers))
    cmd.create(supers,pdbs[0],1,1)
    for idx, (a1, a1resid, a2, a2resid, score) in enumerate(results):
        print(f"\nAlignment {idx+1} (Score: {score})")
        print('Structure 1 - Residues %s - %s' % (a1resid[0],a1resid[-1]))
        print(''.join(a1))
        print('Structure 2 - Residues %s - %s' % (a2resid[0],a2resid[-1]))
        print(''.join(a2))
        cmd.create('mvo',pdbs[0],1,1)
        cmd.super('mvo and resid '+str(a1resid[0])+'-'+str(a1resid[-1]),str(pdbs[1])+' and resid '+str(a2resid[0])+'-'+str(a2resid[-1]))
        cmd.create(supers,'mvo',1,idx+1)
        cmd.delete('mvo')
        
