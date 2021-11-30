import nanome
import tempfile
from nanome.api.structure import Complex
from nanome.util import async_callback, enums, Logs, Process, Color
from nanome.api.shapes import Anchor, Label, Shape
from nanome._internal._structure._bond import _Bond

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
    _atom_group = {}

    @classmethod
    def _fill_atom_valence_table(cls):
        cls._group_to_valence = {1:1, 2:2, 13:3, 14:4, 15:5, 16:6, 17:7, 18:8}
        cls._atom_group["xx"] = 18
        cls._atom_group["h"] = 1
        cls._atom_group["he"] = 18
        cls._atom_group["li"] = 1
        cls._atom_group["be"] = 2
        cls._atom_group["b"] = 13
        cls._atom_group["c"] = 14
        cls._atom_group["n"] = 15
        cls._atom_group["o"] = 16
        cls._atom_group["f"] = 17
        cls._atom_group["ne"] = 18
        cls._atom_group["na"] = 1
        cls._atom_group["mg"] = 2
        cls._atom_group["al"] = 13
        cls._atom_group["si"] = 14
        cls._atom_group["p"] = 15
        cls._atom_group["s"] = 16
        cls._atom_group["cl"] = 17
        cls._atom_group["ar"] = 18
        cls._atom_group["k"] = 1
        cls._atom_group["ca"] = 2
        cls._atom_group["sc"] = 3
        cls._atom_group["ti"] = 4
        cls._atom_group["v"] = 5
        cls._atom_group["cr"] = 6
        cls._atom_group["mn"] = 7
        cls._atom_group["fe"] = 8
        cls._atom_group["co"] = 9
        cls._atom_group["ni"] = 10
        cls._atom_group["cu"] = 11
        cls._atom_group["zn"] = 12
        cls._atom_group["ga"] = 13
        cls._atom_group["ge"] = 14
        cls._atom_group["as"] = 15
        cls._atom_group["se"] = 16
        cls._atom_group["br"] = 17
        cls._atom_group["kr"] = 18
        cls._atom_group["rb"] = 1
        cls._atom_group["sr"] = 2
        cls._atom_group["y"] = 3
        cls._atom_group["zr"] = 4
        cls._atom_group["nb"] = 5
        cls._atom_group["mo"] = 6
        cls._atom_group["tc"] = 7
        cls._atom_group["ru"] = 8
        cls._atom_group["rh"] = 9
        cls._atom_group["pd"] = 10
        cls._atom_group["ag"] = 11
        cls._atom_group["cd"] = 12
        cls._atom_group["in"] = 13
        cls._atom_group["sn"] = 14
        cls._atom_group["sb"] = 15
        cls._atom_group["te"] = 16
        cls._atom_group["i"] = 17
        cls._atom_group["xe"] = 18
        cls._atom_group["cs"] = 1
        cls._atom_group["ba"] = 2
        cls._atom_group["la"] = 3
        cls._atom_group["ce"] = 4
        cls._atom_group["pr"] = 5
        cls._atom_group["nd"] = 6
        cls._atom_group["pm"] = 7
        cls._atom_group["sm"] = 8
        cls._atom_group["eu"] = 9
        cls._atom_group["gd"] = 10
        cls._atom_group["tb"] = 11
        cls._atom_group["dy"] = 12
        cls._atom_group["ho"] = 13
        cls._atom_group["er"] = 14
        cls._atom_group["tm"] = 15
        cls._atom_group["yb"] = 16
        cls._atom_group["lu"] = 17
        cls._atom_group["hf"] = 4
        cls._atom_group["ta"] = 5
        cls._atom_group["w"] = 6
        cls._atom_group["re"] = 7
        cls._atom_group["os"] = 8
        cls._atom_group["ir"] = 9
        cls._atom_group["pt"] = 10
        cls._atom_group["au"] = 11
        cls._atom_group["hg"] = 12
        cls._atom_group["tl"] = 13
        cls._atom_group["pb"] = 14
        cls._atom_group["bi"] = 15
        cls._atom_group["po"] = 16
        cls._atom_group["at"] = 17
        cls._atom_group["rn"] = 18
        cls._atom_group["fr"] = 1
        cls._atom_group["ra"] = 2
        cls._atom_group["ac"] = 3
        cls._atom_group["th"] = 4
        cls._atom_group["pa"] = 5
        cls._atom_group["u"] = 6
        cls._atom_group["np"] = 7
        cls._atom_group["pu"] = 8
        cls._atom_group["am"] = 9
        cls._atom_group["cm"] = 10
        cls._atom_group["bk"] = 11
        cls._atom_group["cf"] = 12
        cls._atom_group["es"] = 13
        cls._atom_group["fm"] = 14
        cls._atom_group["md"] = 15
        cls._atom_group["no"] = 16
        cls._atom_group["lr"] = 17
        cls._atom_group["rf"] = 4
        cls._atom_group["db"] = 5
        cls._atom_group["sg"] = 6
        cls._atom_group["bh"] = 7
        cls._atom_group["hs"] = 8
        cls._atom_group["mt"] = 9
        cls._atom_group["ds"] = 10
        cls._atom_group["rg"] = 11
        cls._atom_group["cn"] = 12
        cls._atom_group["nh"] = 13
        cls._atom_group["fl"] = 14
        cls._atom_group["mc"] = 15
        cls._atom_group["lv"] = 16
        cls._atom_group["ts"] = 17
        cls._atom_group["og"] = 18

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

        self.formal_charges_labels = True

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
        label.font_size = 0.1
        label.color = Color.Black()
        return label

    def bond_type_label(self, bond):
        label = Label()
        anchor = Anchor()
        anchor.anchor_type = nanome.util.enums.ShapeAnchorType.Atom
        anchor.target = bond.atom1.index
        offset = (bond.atom2.position - bond.atom1.position) * 0.5
        anchor.viewer_offset = nanome.util.Vector3(0, 0, -.01)
        anchor.local_offset = offset#nanome.util.Vector3(0, 0, 0)

        label.anchors = [anchor]
        label.text = str(bond.kind).replace("Kind.Covalent", "")
        label.font_size = 0.05
        label.color = Color.Black()
        return label

    @async_callback
    async def add_formal_charge_labels(self, indices):
        deep = await self.request_complexes(indices)
        formal_charges = self.compute_formal_charges(deep)
        labels = []
        for c in deep:
            # for a in c.atoms:
            #     if a.polar_hydrogen:
            #         labels.append(self.formal_charge_label(a, "Polar!"))
            #         a.labelled = True
            # for b in c.bonds:
            #     label = self.bond_type_label(b)
            #     if label.text != "Single":
            #         labels.append(label)
            for a in c.atoms:
                if a in formal_charges:
                    if a.symbol == "H":
                        s_formal_charge = str(formal_charges[a])
                        labels.append(self.formal_charge_label(a, s_formal_charge))
                        a.labelled = True


        Shape.upload_multiple(labels)

    def compute_formal_charges(self, complexes):
        if len(Hydrogens._atom_group) < 1:
            Hydrogens._fill_atom_valence_table()

        formal_charges = {}
        for c in complexes:
            for a in c.atoms:
                atom_symbol = a.symbol.lower()
                if atom_symbol not in Hydrogens._atom_group:
                    Logs.warning("Unkown atom",a.symbol)
                    continue
                group = Hydrogens._atom_group[atom_symbol]
                if group not in Hydrogens._group_to_valence:
                    Logs.warning("Atom",a.symbol,"unknown valence electron")
                    continue
                v_e = Hydrogens._group_to_valence[group]
                n_bonds = self.count_bonded_electrons(a.bonds)
                formal_charges[a] = v_e - n_bonds
        return formal_charges


    def count_bonded_electrons(self, bonds):
        count = 0
        for b in bonds:
            if b.kind == _Bond.Kind.CovalentSingle:
                count += 1
            if b.kind == _Bond.Kind.CovalentDouble:
                count +=2
            if b.kind == _Bond.Kind.CovalentTriple:
                count +=3
        return count
            
                
                
                

def main():
    integrations = [enums.Integrations.hydrogen]
    plugin = nanome.Plugin('Hydrogens', 'A Nanome Plugin to add/remove hydrogens to/from structures', 'Hydrogens', True, integrations=integrations)
    plugin.set_plugin_class(Hydrogens)
    plugin.run()

if __name__ == '__main__':
    main()
