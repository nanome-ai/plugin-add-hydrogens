import nanome
import tempfile
from nanome.api.structure import Complex
from nanome.util import async_callback, enums, Logs, Process

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
        self.input_file = tempfile.NamedTemporaryFile(delete=False, suffix='.sdf', dir=self.temp_dir.name)
        self.output_file = tempfile.NamedTemporaryFile(delete=False, suffix='.sdf', dir=self.temp_dir.name)
        self.complex_names = []

        self.set_plugin_list_button(enums.PluginListButtonType.advanced_settings, 'pH Settings')
        self.integration.hydrogen_add = self.integration_add
        self.integration.hydrogen_remove = self.integration_remove
        self.ph = '7.4'

        self.create_settings_menu()

    def create_settings_menu(self):
        menu = nanome.ui.Menu()
        self.menu = menu

        menu.title = 'Settings'
        menu.width = 0.5
        menu.height = 0.2

        menu.root.layout_orientation = menu.root.LayoutTypes.horizontal
        menu.root.padding_type = menu.root.PaddingTypes.ratio
        menu.root.set_padding(top=0.25, down=0.25, left=0.05, right=0.05)

        ln_lbl = menu.root.create_child_node()
        lbl = ln_lbl.add_new_label('pH')
        lbl.text_horizontal_align = enums.HorizAlignOptions.Middle
        lbl.text_vertical_align = enums.VertAlignOptions.Middle

        ln_inp = menu.root.create_child_node()
        ln_inp.forward_dist = 0.001
        inp = ln_inp.add_new_text_input()
        inp.number = True
        inp.input_text = self.ph

        def change_ph(input):
            self.ph = input.input_text
        inp.register_changed_callback(change_ph)

    def on_advanced_settings(self):
        self.menu.enabled = True
        self.update_menu(self.menu)

    @async_callback
    async def on_run(self):
        shallow = await self.request_complex_list()
        indices_selected = [c.index for c in shallow if c.get_selected()]

        if not indices_selected:
            self.send_notification(enums.NotificationTypes.warning, 'Please select a complex.')
            return

        self.set_plugin_list_button(enums.PluginListButtonType.run, 'Running...', False)

        # get selected complexes and add hydrogens
        complexes = await self.request_complexes(indices_selected)
        await self.add_hydrogens(complexes)
        self.update_structures_deep(complexes)

        self.set_plugin_list_button(enums.PluginListButtonType.run, 'Run', True)
        self.send_notification(enums.NotificationTypes.success, f'Hydrogens calculated with pH {self.ph}')

    def on_stop(self):
        self.temp_dir.cleanup()

    @async_callback
    async def integration_add(self, request):
        complexes = request.get_args()
        await self.add_hydrogens(complexes)
        request.send_response(complexes)

    def integration_remove(self, request):
        complexes = request.get_args()
        self.remove_hydrogens(complexes)
        request.send_response(complexes)

    async def add_hydrogens(self, complexes):
        Logs.debug('Add H')

        # remove all hydrogens before beginning
        self.remove_hydrogens(complexes)
        self.complex_names = []

        for complex in complexes:
            self.complex_names.append(complex.name)
            complex.io.to_sdf(self.input_file.name)

            # remember atoms by position
            atom_by_position = dict()
            for atom in complex.atoms:
                atom_by_position[get_position_key(atom)] = atom

            # compute all hydrogens
            result_complex = await self.compute_hydrogens()
            if not result_complex:
                continue

            # add hydrogens to original complex
            self.match_and_update(atom_by_position, result_complex)

            # compute polar hydrogens
            result_complex = await self.compute_hydrogens(True)
            if not result_complex:
                continue

            # mark hydrogens as polar in original complex
            self.match_and_update(atom_by_position, result_complex, True)

    def remove_hydrogens(self, complexes):
        Logs.debug('Remove H')

        for complex in complexes:
            for atom in list(complex.atoms):
                if atom.symbol != 'H' or not atom.selected:
                    continue

                residue = atom.residue
                residue.remove_atom(atom)
                for bond in list(atom.bonds):
                    residue.remove_bond(bond)

    async def compute_hydrogens(self, polar=False):
        p = Process()
        p.executable_path = 'nanobabel'
        p.args = ['hydrogen', '-add', '-ph', self.ph, '-i', self.input_file.name, '-o', self.output_file.name]
        p.output_text = True
        p.on_error = Logs.warning
        p.on_output = Logs.debug

        if polar:
            p.args += ['-pl']

        # if error, notify user and return
        exit_code = await p.start()
        if exit_code:
            self.send_notification(enums.NotificationTypes.error, 'Error computing hydrogens')
            return

        return Complex.io.from_sdf(path=self.output_file.name)

    def match_and_update(self, atom_by_position, result_complex, polar=False):
        """
        Match and add hydrogens to original complex using positions to match.

        :param atom_by_position: dict mapping position key to atom in source complex
        :type atom_by_position: dict
        :param result_complex: Output complex from hydrogens calculation
        :type result_complex: :class:`nanome.structure.Complex`
        :param polar: If true, update hydrogens as polar, otherwise add hydrogens to source complex, defaults to False
        :type polar: bool, optional
        """

        num_non_hydrogens = 0
        num_unknown_atoms = 0

        for atom in result_complex.atoms:
            if atom.symbol != 'H':
                num_non_hydrogens += 1
                continue

            # H can only have 1 bond
            bond = next(atom.bonds)
            bonded_atom = bond.atom1 if bond.atom1.symbol != 'H' else bond.atom2
            bonded_atom_key = get_position_key(bonded_atom)

            if bonded_atom_key not in atom_by_position:
                num_unknown_atoms += 1
                continue

            source_atom = atom_by_position[bonded_atom_key]
            if not source_atom.selected:
                continue

            if not polar:
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

            else:
                for bond in source_atom.bonds:
                    if bond.atom1.symbol == 'H':
                        bond.atom1.polar_hydrogen = True
                    if bond.atom2.symbol == 'H':
                        bond.atom2.polar_hydrogen = True

        if num_unknown_atoms:
            Logs.warning(f'{self.complex_names}: {num_unknown_atoms}/{num_non_hydrogens} atoms not found')


def main():
    integrations = [enums.Integrations.hydrogen]
    plugin = nanome.Plugin('Hydrogens', 'A Nanome Plugin to add/remove hydrogens to/from structures', 'Hydrogens', True, integrations=integrations)
    plugin.set_plugin_class(Hydrogens)
    plugin.run()

if __name__ == '__main__':
    main()
