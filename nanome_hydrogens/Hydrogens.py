from .ComplexUtils import ComplexUtils as utils
import nanome
from nanome.util import async_callback, Logs, enums
from nanome.util.enums import Integrations, NotificationTypes
from rdkit import Chem
import tempfile

def get_position_key(atom):
    """
    Get atom position as tuple of integers to be used as lookup key.
    Rounds the coordinates to 4 decimal places before multiplying by 50
    to get unique integer-space coordinates, and avoid floating point errors.
    :param atom: Atom to get key for
    :type atom: :class:`nanome.structure.Atom`
    :return: Position key tuple
    :rtype: (int, int, int)
    """
    return tuple(map(lambda x: int(50 * round(x, 4)), atom.position))

class Hydrogens(nanome.AsyncPluginInstance):
    def start(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.input_file_rm = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf", dir=self.temp_dir.name)
        self.output_file_rm = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf", dir=self.temp_dir.name)
        self.input_file_add = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf", dir=self.temp_dir.name)
        self.output_file_add = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf", dir=self.temp_dir.name)
        self.processes = []
        self.integration.hydrogen_add = self.add_hydrogens
        self.integration.hydrogen_remove = self.remove_hydrogens

    @async_callback
    async def on_run(self):
        shallow = await self.request_complex_list()
        indices_selected = [c.index for c in shallow if c.get_selected()]

        if not indices_selected:
            self.send_notification(enums.NotificationTypes.warning, 'Please select a complex.')
            return
        
        self.set_plugin_list_button(enums.PluginListButtonType.run, 'Running...', False)

        deep = await self.request_complexes(indices_selected)
        result = await self.add_hydrogens(complexes=deep)
        self.update_structures_deep(result)

        self.set_plugin_list_button(enums.PluginListButtonType.run, 'Run', True)
        self.send_notification(enums.NotificationTypes.success, f'Hydrogens calculated')


    @async_callback
    async def add_hydrogens(self, request=None, complexes=None):
        Logs.debug('Add Hs')

        if request:
            complexes = request.get_args()
        
        for complex in complexes:
            # remember atoms by position
            atom_by_position = dict()
            for atom in complex.atoms:
                atom_by_position[get_position_key(atom)] = atom

            # nanome complex -> sdf -> rdkit molecule
            complex.io.to_sdf(self.input_file_add.name)
            rdmol = Chem.SDMolSupplier(self.input_file_add.name, removeHs=True)[0]
            # +Hs
            rdmol = Chem.AddHs(rdmol, addCoords=True, addResidueInfo=False)
            # rdkit molecule -> sdf -> nanome complex
            sdw = Chem.SDWriter(self.output_file_add.name)
            sdw.write(rdmol)
            sdw.close()
            H_complex = nanome.structure.Complex.io.from_sdf(path=self.output_file_add.name)
            
            if not H_complex:
                continue
            # add hydrogens to original complex
            self.match_and_update(atom_by_position, H_complex)

        if request:
            request.send_response(complexes)

        return complexes

    def remove_hydrogens(self, request=None, complexes=None):
        Logs.debug('Remove H')

        if request:
            complexes = request.get_args()

        for complex in complexes:
            for atom in list(complex.atoms):
                if atom.symbol != 'H' or not atom.selected:
                    continue

                residue = atom.residue
                residue.remove_atom(atom)
                for bond in list(atom.bonds):
                    residue.remove_bond(bond)

        if request:
            request.send_response(complexes)

        return complexes

    def match_and_update(self, atom_by_position, result_complex):
        """
        Match and add hydrogens to original complex using positions to match.
        :param atom_by_position: dict mapping position key to atom in source complex
        :type atom_by_position: dict
        :param result_complex: Output complex from hydrogens calculation
        :type result_complex: :class:`nanome.structure.Complex`
        """
        for atom in result_complex.atoms:
            if atom.symbol != 'H':
                continue

            # H can only have 1 bond
            bond = next(atom.bonds)
            bonded_atom = bond.atom1 if bond.atom1.symbol != 'H' else bond.atom2
            bonded_atom_key = get_position_key(bonded_atom)

            if bonded_atom_key not in atom_by_position:
                Logs.warning(f'H {atom.serial} bonded with unknown atom {bonded_atom.symbol} at {bonded_atom.position}')
                continue

            source_atom = atom_by_position[bonded_atom_key]
            if not source_atom.selected:
                continue

            new_atom = atom._shallow_copy()
            new_atom._display_mode = source_atom._display_mode
            new_atom.is_het = source_atom.is_het
            new_atom.selected = source_atom.selected

            # copy bfactor and occupancy for surface coloring
            new_atom.bfactor = source_atom.bfactor
            new_atom.occupancy = source_atom.occupancy

            new_bond = bond._shallow_copy()
            new_bond.atom1 = source_atom
            new_bond.atom2 = new_atom

            residue = source_atom.residue
            residue.add_atom(new_atom)
            residue.add_bond(new_bond)

    def on_stop(self):
        self.temp_dir.cleanup()

def main():
    plugin = nanome.Plugin('Hydrogens', 'A nanome integration plugin to add and remove hydrogens to/from structures', 'Hydrogens', False, integrations=[Integrations.hydrogen])
    plugin.set_plugin_class(Hydrogens)
    plugin.run()


if __name__ == '__main__':
    main()