from pathlib import Path

from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.Structure import Structure


def read_pdb(path: Path, sid: str = "structure") -> Structure:
    parser = PDBParser(QUIET=True)
    return parser.get_structure(sid, str(path))


def get_chain(structure: Structure, chain_id: str) -> Chain:
    """Note: Ligands, water molecules, ions, and other molecules are not removed."""
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                return chain
    raise ValueError(f"No chain with id={chain_id} found.")


def save_chain(path: str, chain: Chain, name: str = "Structure"):
    new_chain = Chain(chain.id)
    for res in chain:
        if is_aa(res):
            new_chain.add(res)
    structure = Structure(name)
    model = Model(0)
    structure.add(model)
    model.add(new_chain)
    io = PDBIO()
    io.set_structure(structure)
    io.save(path)
