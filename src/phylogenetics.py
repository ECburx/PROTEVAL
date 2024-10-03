import os
import re
import subprocess
from pathlib import Path
from typing import List

import dendropy
import numpy as np
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from joblib import Parallel, delayed
from tqdm import tqdm


def qscore(pdb1, pdb2, r0: int = 4) -> float:
    """
    Q_score  = Q_length * Q_shape
    Q_length = N_align^2 / (N1 * N2)
    Q_shape  = 1 / (1 + (RMSD / R0)^2)
    https://doi.org/10.1107/s0907444904026460
    https://doi-org.ezp.lib.cam.ac.uk/10.1093/molbev/msaa100
    :param r0: an empirical parameter (chosen at 4Å)
    :return:
    """
    tm = tmalign(pdb1, pdb2)
    q_length = (tm["Nalign"] ** 2) / (tm["N1"] * tm["N2"])
    q_shape = 1 / (1 + (tm["rmsd"] / r0) ** 2)
    return q_length * q_shape


def tmalign(pdb1_path: str, pdb2_path: str, tmalign_path: str = "../tmalign/TMalign") -> dict:
    """
    Wrapper of TMalign.
    :return:
        {'N1':              Length of chain 1,
         'N2':              Length of chain 2,
         'Nalign':          Number of aligned residues,
         'rmsd':            RMSD,
         'mean_TMscore':    Mean TM-score:
                            (score normalized by length of Chain_1 + score normalized by length of Chain_2) / 2
        }
    """
    if not os.path.exists(pdb1_path):
        raise FileNotFoundError(pdb1_path)
    if not os.path.exists(pdb2_path):
        raise FileNotFoundError(pdb2_path)

    result = {}
    process = subprocess.run(
        f"{tmalign_path} {pdb1_path} {pdb2_path}",
        shell=True,
        capture_output=True,
        text=True
    )
    output = process.stdout

    match_len_c1 = re.search(r'Length of Chain_1:\s*(\d+)', output)
    match_len_c2 = re.search(r'Length of Chain_2:\s*(\d+)', output)
    match_Nalign = re.search(r'Aligned length=\s*(\d+)', output)
    match_tmRMSD = re.search(r'RMSD=\s*([\d\.]+)', output)
    match_score1 = re.search(r'TM-score=\s*([\d\.]+) \(if normalized by length of Chain_1', output)
    match_score2 = re.search(r'TM-score=\s*([\d\.]+) \(if normalized by length of Chain_2', output)

    if match_len_c1: result['N1'] = int(match_len_c1.group(1))
    if match_len_c2: result['N2'] = int(match_len_c2.group(1))
    if match_Nalign: result['Nalign'] = int(match_Nalign.group(1))
    if match_tmRMSD: result['rmsd'] = float(match_tmRMSD.group(1))
    if match_score1 and match_score2:
        result['mean_TMscore'] = (float(match_score1.group(1)) + float(match_score2.group(1))) / 2

    if len(result.keys()) != 5:
        raise Warning(f"TM-align might not work as expected.\n"
                      f"{tmalign_path} {pdb1_path} {pdb2_path}\n"
                      f"Check the following output:\n"
                      f"stdout:\n"
                      f"{output}\n"
                      f"stderr:\n"
                      f"{process.stderr}")
    return result


def rmsd_matrix(pdb_list: List[str]):
    n_pdb = len(pdb_list)
    matrix = np.zeros((n_pdb, n_pdb))
    for i in tqdm(range(n_pdb)):
        for j in range(n_pdb):
            if i == j:
                matrix[i, j] = 0
            else:
                matrix[i, j] = tmalign(pdb_list[i], pdb_list[j])["rmsd"]
    return matrix


def qs_phylogenetic_tree(pdb_list: List[str], names: List[str], n_jobs: int = 10):
    """
    Generate a phylogenetic tree from Q-scores of a list of PDB files.
    :param pdb_list: A list of PDB file paths.
    :param names: A list of names (identifiers) of the PDB files.
    :param n_jobs
    :return: {"dist_matrix", "qscore_matrix", "tree", "newick"}
    """
    # Q-score Matrix.
    n_pdb = len(pdb_list)
    qscore_matrix = np.zeros((n_pdb, n_pdb))

    # for i in tqdm(range(n_pdb)):
    #     for j in range(i, n_pdb):
    #         if i == j:
    #             qscore_matrix[i, j] = 1.0
    #         else:
    #             score_ij = qscore(pdb_list[i], pdb_list[j])
    #             score_ji = qscore(pdb_list[j], pdb_list[i])
    #             score_avg = (score_ij + score_ji) / 2.0
    #             qscore_matrix[i, j] = score_avg
    #             qscore_matrix[j, i] = score_avg

    def compute_score(i, j):
        if i == j:
            return i, j, 1.0
        else:
            score_ij = qscore(pdb_list[i], pdb_list[j])
            score_ji = qscore(pdb_list[j], pdb_list[i])
            score_avg = (score_ij + score_ji) / 2.0
            return i, j, score_avg

    results = Parallel(n_jobs=n_jobs)(
        delayed(compute_score)(i, j)
        for i, j in tqdm([
            (i, j)
            for i in range(n_pdb)
            for j in range(i, n_pdb)
        ])
    )

    for i, j, score in results:
        qscore_matrix[i, j] = score
        qscore_matrix[j, i] = score

    # Distance Matrix.
    dist_matrix = 1 - qscore_matrix

    # Phylo Tree
    phylo_tree = DistanceTreeConstructor().upgma(
        distance_matrix=DistanceMatrix(
            names=names,
            matrix=[dist_matrix[i, :i + 1].tolist() for i in range(len(dist_matrix))]
        )
    )

    return {
        "dist_matrix": dist_matrix,
        "qscore_matrix": qscore_matrix,
        "tree": phylo_tree,
        "newick": phylo_tree.format("newick")
    }


def normalize(input_path: Path, output_path: Path, schema: str = "newick", mode: str = "sum"):
    def normalize_tree(t):
        assert mode in ["sum", "max"]
        branch_lengths = [edge.length for edge in t.edges() if edge.length is not None]
        if not branch_lengths:
            return t

        sum_len = sum(branch_lengths)
        max_len = max(branch_lengths)
        for edge in t.edges():
            if edge.length is None:
                continue
            if mode == "sum":
                edge.length /= sum_len
            if mode == "max":
                edge.length /= max_len
        return t

    input_trees = dendropy.TreeList.get(path=input_path, schema=schema)
    normalized_trees = dendropy.TreeList()

    for tree in input_trees:
        normalized_trees.append(normalize_tree(tree))

    normalized_trees.write(path=output_path, schema=schema)
