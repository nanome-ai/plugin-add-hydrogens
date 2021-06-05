from os.path import dirname, join
import nanome

BASE_DIR = dirname(__file__)
ELECTRONEGATIVITIES = {}

class ComplexUtils:

    with open(join(BASE_DIR, 'electronegativities.csv'), 'r') as file:
        lines = file.readlines()
        for line in lines:
            element, eN = line.split(',')
            ELECTRONEGATIVITIES[element] = float(eN)

    @staticmethod
    def align_to(complex, reference_complex):
        m = complex.get_complex_to_workspace_matrix()
        for atom in complex.atoms:
            atom.position = m * atom.position
        complex.old_position = complex.position
        complex.old_rotation = complex.rotation
        complex.position = reference_complex.position
        complex.rotation = reference_complex.rotation
        m = complex.get_workspace_to_complex_matrix()
        for atom in complex.atoms:
            atom.position = m * atom.position

    @staticmethod
    def combine_ligands(receptor, ligands):
        combined_ligands = nanome.structure.Complex()
        combined_ligands.names = []
        for ligand in ligands:
            ComplexUtils.align_to(ligand, receptor)
            combined_ligands.names.append(ligand.full_name)
            for molecule in ligand.molecules:
                combined_ligands.add_molecule(molecule)
        return combined_ligands

    @staticmethod
    def convert_to_conformers(complexes):
        for i in range(len(complexes)):
            complex_index = complexes[i].index
            complexes[i] = complexes[i].convert_to_conformers()
            complexes[i].index = complex_index

    @staticmethod
    def convert_to_frames(complexes):
        for i in range(len(complexes)):
            complex_index = complexes[i].index
            complexes[i] = complexes[i].convert_to_frames()
            complexes[i].index = complex_index

    @staticmethod
    def markPolarHydrogens(complex):
        for atom in complex.atoms:
            if atom.polar_hydrogen: continue
            for bond in atom.bonds:
                atom1, atom2 = bond.atom1, bond.atom2
                eN1 = ELECTRONEGATIVITIES[atom1.symbol]
                eN2 = ELECTRONEGATIVITIES[atom2.symbol]
                eN_diff = abs(eN1-eN2)
                if (eN_diff > 0.4 and eN_diff <= 1.8):
                    atom1.polar_hydrogen = atom1.symbol == 'H'
                    atom2.polar_hydrogen = atom1.symbol == 'H'

    @staticmethod
    def reidentify(target, source):
        target.index = source.index
        target.name = source.name
        target.position = source.position
        target.rotation = source.rotation

    @staticmethod
    def reset_transform(complex):
        complex.position = complex.old_position
        complex.rotation = complex.old_rotation