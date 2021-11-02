import nanome
import tempfile
from nanome.api.structure import Complex
from nanome.util import async_callback, enums, Logs, Process, Color
from nanome.api.shapes import Anchor, Label, Shape

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

        self.set_plugin_list_button(enums.PluginListButtonType.advanced_settings, 'pH Settings')
        self.integration.hydrogen_add = self.add_hydrogens
        self.integration.hydrogen_remove = self.remove_hydrogens
        self.ph = '7.4'

        self.create_settings_menu()

    def create_settings_menu(self):
        menu = nanome.ui.Menu()
        self.menu = menu

        menu.title = 'Settings'
        menu.width = 0.8
        menu.height = 0.8

        menu.root.layout_orientation = menu.root.LayoutTypes.vertical
        menu.root.padding_type = menu.root.PaddingTypes.ratio
        menu.root.set_padding(top=0.05, down=0.05, left=0.05, right=0.05)

        main = menu.root.create_child_node()
        main.layout_orientation = main.LayoutTypes.horizontal

        ln_lbl = main.create_child_node()
        lbl = ln_lbl.add_new_label('pH')
        lbl.text_horizontal_align = enums.HorizAlignOptions.Left
        lbl.text_vertical_align = enums.VertAlignOptions.Middle

        ln_inp = main.create_child_node()
        ln_inp.forward_dist = 0.001
        inp = ln_inp.add_new_text_input()
        inp.number = True
        inp.input_text = self.ph

        main2 = menu.root.create_child_node()
        main2.layout_orientation = main2.LayoutTypes.horizontal

        ln_lbl2 = main2.create_child_node()
        ln_lbl2.text_horizontal_align = enums.HorizAlignOptions.Left
        ln_lbl2.text_vertical_align = enums.VertAlignOptions.Middle
        lbl2 = ln_lbl2.add_new_label("Formal\ncharges\nlabels")
        lbl2.text_auto_size = False
        lbl2.text_size = 0.5

        ln_lbl3 = main2.create_child_node()
        ln_lbl3.padding_type = ln_lbl3.PaddingTypes.fixed
        ln_lbl3.padding = (0.1, 0.1, 0.1, 0.1)
        lbl3 = ln_lbl3.add_new_toggle_switch("")
        lbl3.forward_dist = 0.001
        lbl3.selected = False

        self.formal_charges_labels = False

        def change_ph(input):
            self.ph = input.input_text
        inp.register_changed_callback(change_ph)

        def set_formal_charges(input):
            self.formal_charges_labels = input.selected

        lbl3.register_pressed_callback(set_formal_charges)

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
        deep = await self.request_complexes(indices_selected)
        result = await self.add_hydrogens(complexes=deep)
        await self.update_structures_deep(result)

        if self.formal_charges_labels:
            self.add_formal_charge_labels(indices_selected)

        self.set_plugin_list_button(enums.PluginListButtonType.run, 'Run', True)
        self.send_notification(enums.NotificationTypes.success, f'Hydrogens calculated with pH {self.ph}')

    def on_stop(self):
        self.temp_dir.cleanup()

    @async_callback
    async def add_hydrogens(self, request=None, complexes=None):
        Logs.debug('Add H')

        if request:
            complexes = request.get_args()

        # remove all hydrogens before beginning
        self.remove_hydrogens(complexes=complexes)

        for complex in complexes:
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

    async def compute_hydrogens(self, polar=False):
        p = Process()
        p.executable_path = 'nanobabel'
        p.args = ['hydrogen', '-add', '-ph', self.ph, '-i', self.input_file.name, '-o', self.output_file.name]
        p.output_text = True
        p.on_error = Logs.error
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

    def formal_charge_label(self, atom, text):
        label = Label()
        anchor = Anchor()
        anchor.anchor_type = nanome.util.enums.ShapeAnchorType.Atom
        anchor.target = atom.index
        anchor.viewer_offset = nanome.util.Vector3(0, 0, -.01)
        anchor.local_offset = nanome.util.Vector3(0, 0, 0)

        label.anchors = [anchor]
        label.text = text
        label.font_size = 0.05
        label.color = Color.Black()
        return label

    @async_callback
    async def add_formal_charge_labels(self, indices):
        deep = await self.request_complexes(indices)
        labels = []
        for c in deep:
            for a in c.atoms:
                if a.polar_hydrogen:
                    labels.append(self.formal_charge_label(a, "Polar!"))
        Shape.upload_multiple(labels)

def main():
    integrations = [enums.Integrations.hydrogen]
    plugin = nanome.Plugin('Hydrogens', 'A Nanome Plugin to add/remove hydrogens to/from structures', 'Hydrogens', True, integrations=integrations)
    plugin.set_plugin_class(Hydrogens)
    plugin.run()

if __name__ == '__main__':
    main()
