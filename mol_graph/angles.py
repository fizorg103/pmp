import os
from rdkit import Chem
from mol_graph.xyz2mol import read_xyz_file, xyz2mol
import numpy as np


def get_dihedral_deg(x):
    file_name = 'input'
    file_name = os.path.join(file_name, 'structures')
    random_molecule = x['molecule_name']
    file_name = os.path.join(file_name, random_molecule + '.xyz')
    atomicNumList, charge, xyz_coordinates = read_xyz_file(file_name)

    charged_fragments = True
    quick = True

    mol = xyz2mol(atomicNumList, charge, xyz_coordinates, charged_fragments, quick)

    atom_id_0 = x['atom_index_0']
    atom_id_1 = x['atom_index_1']

    atom_0 = mol.GetAtomWithIdx(atom_id_0)
    atom_0_neig = atom_0.GetNeighbors()[0].GetIdx()

    atom_1 = mol.GetAtomWithIdx(atom_id_1)
    atom_1_neig = atom_1.GetNeighbors()[0].GetIdx()

    conf = mol.GetConformer()

    angle = Chem.rdMolTransforms.GetDihedralDeg(conf, atom_id_0, atom_0_neig, atom_1_neig, atom_id_1)
    return np.abs(angle)
