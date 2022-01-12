import nanome
import tempfile
from nanome.api.structure import Complex
from nanome.util import async_callback, enums, Logs, Process, Color
from nanome.api.shapes import Anchor, Label, Shape
from nanome._internal._structure._bond import _Bond

from rdkit import Chem
from rdkit.Chem import AllChem
import psi4
import resp

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
    _valence_atoms = {}

    @classmethod
    def _fill_atom_valence_table(cls):
        cls._valence_atoms["h"] =  1
        cls._valence_atoms["li"] =  1
        cls._valence_atoms["be"] =  2
        cls._valence_atoms["b"] =  3
        cls._valence_atoms["c"] =  4
        cls._valence_atoms["n"] =  5
        cls._valence_atoms["o"] =  6
        cls._valence_atoms["f"] =  7
        cls._valence_atoms["ne"] =  8
        cls._valence_atoms["na"] =  1
        cls._valence_atoms["mg"] =  2
        cls._valence_atoms["al"] =  3
        cls._valence_atoms["si"] =  4
        cls._valence_atoms["p"] =  5
        cls._valence_atoms["s"] =  6
        cls._valence_atoms["cl"] =  7
        cls._valence_atoms["ar"] =  8
        cls._valence_atoms["k"] =  1
        cls._valence_atoms["ca"] =  2
        cls._valence_atoms["ga"] =  3
        cls._valence_atoms["ge"] =  4
        cls._valence_atoms["as"] =  5
        cls._valence_atoms["se"] =  6
        cls._valence_atoms["br"] =  7
        cls._valence_atoms["kr"] =  8
        cls._valence_atoms["rb"] =  1
        cls._valence_atoms["sr"] =  2
        cls._valence_atoms["in"] =  3
        cls._valence_atoms["sn"] =  4
        cls._valence_atoms["sb"] =  5
        cls._valence_atoms["te"] =  6
        cls._valence_atoms["i"] =  7
        cls._valence_atoms["xe"] =  8
        cls._valence_atoms["cs"] =  1
        cls._valence_atoms["ba"] =  2
        cls._valence_atoms["tl"] =  3
        cls._valence_atoms["pb"] =  4
        cls._valence_atoms["bi"] =  5
        cls._valence_atoms["po"] =  6
        cls._valence_atoms["at"] =  7
        cls._valence_atoms["rn"] =  8
        cls._valence_atoms["fr"] =  1
        cls._valence_atoms["ra"] =  2
        cls._valence_atoms["nh"] =  3
        cls._valence_atoms["fl"] =  4
        cls._valence_atoms["mc"] =  5
        cls._valence_atoms["lv"] =  6
        cls._valence_atoms["ts"] =  7
        cls._valence_atoms["og"] =  8
        cls._valence_atoms["he"] =  2

    def start(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.input_file = tempfile.NamedTemporaryFile(delete=False, suffix='.sdf', dir=self.temp_dir.name)
        self.output_file = tempfile.NamedTemporaryFile(delete=False, suffix='.sdf', dir=self.temp_dir.name)

        self.set_plugin_list_button(enums.PluginListButtonType.advanced_settings, 'pH Settings')
        self.integration.hydrogen_add = self.add_hydrogens
        self.integration.hydrogen_remove = self.remove_hydrogens
        self.ph = '7.4'

        self.formal_labels = []
        self.partial_labels = []

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
        inp.sizing_value = 0.5

        main2 = menu.root.create_child_node()
        main2.layout_orientation = main2.LayoutTypes.horizontal

        ln_lbl2 = main2.create_child_node()
        ln_lbl2.text_horizontal_align = enums.HorizAlignOptions.Middle
        ln_lbl2.text_vertical_align = enums.VertAlignOptions.Middle
        ln_lbl2.padding_type = ln_lbl2.PaddingTypes.fixed
        ln_lbl2.padding = (0.0, 0.0, 0.2, 0.0)
        lbl2 = ln_lbl2.add_new_label("Formal charges labels")
        lbl2.text_auto_size = False
        lbl2.text_size = 0.3

        ln_lbl4 = main2.create_child_node()
        ln_lbl4.padding_type = ln_lbl4.PaddingTypes.fixed
        ln_lbl4.padding = (0.0, 0.0, 0.2, 0.0)
        ln_lbl4.text_horizontal_align = enums.HorizAlignOptions.Middle
        ln_lbl4.text_vertical_align = enums.VertAlignOptions.Middle
        lbl4 = ln_lbl4.add_new_label("Partial charges labels")
        lbl4.text_auto_size = False
        lbl4.text_size = 0.3

        main3 = menu.root.create_child_node()
        main3.layout_orientation = main3.LayoutTypes.horizontal

        ln_lbl3 = main3.create_child_node()
        ln_lbl3.padding_type = ln_lbl3.PaddingTypes.fixed
        ln_lbl3.padding = (0.0, 0.1, 0.1, 0.0)
        lbl3 = ln_lbl3.add_new_toggle_switch("")
        ln_lbl3.forward_dist = 0.001
        lbl3.selected = False

        ln_lbl5 = main3.create_child_node()
        ln_lbl5.padding_type = ln_lbl5.PaddingTypes.fixed
        ln_lbl5.padding = (0.0, 0.1, 0.1, 0.0)
        lbl5 = ln_lbl5.add_new_toggle_switch("")
        ln_lbl5.forward_dist = 0.001
        lbl5.selected = False
        

        self.formal_charges_labels = False
        self.partial_charges_labels = False

        def change_ph(input):
            self.ph = input.input_text
        inp.register_changed_callback(change_ph)

        def set_formal_charges(input):
            if len(self.formal_labels) > 0:
                Shape.destroy_multiple(self.formal_labels)
                self.formal_labels = []
            self.formal_charges_labels = input.selected

        def set_partial_charges(input):
            if len(self.partial_labels) > 0:
                Shape.destroy_multiple(self.partial_labels)
                self.partial_labels = []
            self.partial_charges_labels = input.selected

        lbl3.register_pressed_callback(set_formal_charges)
        lbl5.register_pressed_callback(set_partial_charges)

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
            self.add_formal_charges_labels(indices_selected)
        if self.partial_charges_labels:
            self.add_partial_charges_labels(indices_selected)

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

    def charge_label(self, atom, text):
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
    async def add_formal_charges_labels(self, indices):
        deep = await self.request_complexes(indices)
        formal_charges = self.compute_formal_charges(deep)
        self.formal_labels = []
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
                    if formal_charges[a] != 0:
                        s_formal_charge = str(formal_charges[a])
                        self.formal_labels.append(self.charge_label(a, s_formal_charge))
                        a.labelled = True

        if len(self.formal_labels) > 0:
            self.send_notification(enums.NotificationTypes.message, 'Adding '+str(len(self.formal_labels))+' labels')
            Shape.upload_multiple(self.formal_labels)
        else:
            self.send_notification(enums.NotificationTypes.warning, 'No formal charges label added')


    @async_callback
    async def add_partial_charges_labels(self, indices):
        deep = await self.request_complexes(indices)
        partial_charges = self.compute_partial_charges(deep, "mmff")
        # partial_charges = self.compute_partial_charges(deep, "gasteiger")
        # partial_charges = self.compute_partial_charges(deep, "psi4")
        self.partial_labels = []
        for c in deep:
            idA = 0
            p_c = partial_charges[c.index]
            for a in c.atoms:
                if p_c[idA] != 0:
                    s_partial_charge = str(round(p_c[idA], 3))
                    self.partial_labels.append(self.charge_label(a, s_partial_charge))
                    a.labelled = True
                idA+=1

        if len(self.partial_labels) > 0:
            self.send_notification(enums.NotificationTypes.message, 'Adding '+str(len(self.partial_labels))+' labels')
            Shape.upload_multiple(self.partial_labels)
        else:
            self.send_notification(enums.NotificationTypes.warning, 'No partial charges label added')

    def compute_formal_charges(self, complexes):
        if len(Hydrogens._valence_atoms) < 1:
            Hydrogens._fill_atom_valence_table()

        formal_charges = {}
        for c in complexes:
            for a in c.atoms:
                atom_symbol = a.symbol.lower()
                if atom_symbol not in Hydrogens._valence_atoms:
                    Logs.warning("Unkown atom",a.symbol)
                    continue
                
                va_e = Hydrogens._valence_atoms[atom_symbol]
                n_bonds = self.count_bonded_electrons(a.bonds)
                
                if atom_symbol == "h":
                    non_b_e = va_e - n_bonds
                else:
                    #Not D orbital (most of cases)
                    
                    non_b_e = (8 - va_e) - n_bonds
                    
                    #D orbital (Phosphate / Sulfate / Silica ...)
                    if n_bonds > (8 - va_e) and non_b_e%2 == 0 and atom_symbol in ["p", "s", "ar", "si", "kr", "xe"]:
                        non_b_e = 0

                formal_charges[a] = -non_b_e

        return formal_charges

    def compute_partial_charges(self, complexes, method):
        partial_charges = {}
        Logs.message("Computing partial charges with method '", method+"'")

        for c in complexes:
            molecule = c._molecules[c.current_frame]
            temp_sdf = tempfile.NamedTemporaryFile(delete=False, suffix='.sdf')
            c.io.to_sdf(temp_sdf.name)

            rdmol = next(Chem.SDMolSupplier(temp_sdf.name, removeHs=False, sanitize=False))

            if method == "gasteiger":
                rdmol.UpdatePropertyCache(strict=False)
                Chem.rdmolops.FastFindRings(rdmol)
                AllChem.ComputeGasteigerCharges(rdmol)
                charges = [a.GetDoubleProp('_GasteigerCharge') for a in rdmol.GetAtoms()]
                partial_charges[c.index] = charges
            elif method == "mmff":
                rdmol.UpdatePropertyCache(strict=False)
                Chem.rdmolops.FastFindRings(rdmol)
                mp = AllChem.MMFFGetMoleculeProperties(rdmol)
                charges = [mp.GetMMFFPartialCharge(i) for i in range(rdmol.GetNumAtoms())]
                partial_charges[c.index] = charges

            # TODO: fix this
            # elif method == "psi4" or method == "resp":
            #     if len(list(c.atoms)) > 30:
            #         self.send_notification(enums.NotificationTypes.warning, 'Psi4 partial charges method will take a while...')

            #     xyz = Chem.rdmolfiles.MolToXYZBlock(rdmol)
            #     partial_charges[c.index] = psi4Charges(xyz)


        return partial_charges


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


#From https://github.com/hesther/espsim
def psi4Charges(xyz,
                    basisPsi4 = '3-21G',
                    methodPsi4 = 'scf',
                    gridPsi4 = 1):
    """
    Calculates RESP charges via Psi4.
    :param xyz: String of xyz file of an embedded molecule.
    :param basisPsi4: (optional) Basis set.
    :param methodPsi4: (optional) Method.
    :param gridPsi4: (optional) Integer grid point density for ESP evaluation.
    :return: Array of RESP partial charges.
    """
    Logs.message("Running psi4 resp...")
    
    mol = psi4.core.Molecule.from_string(xyz, dtype='xyz')
    mol.update_geometry()

    options = {'VDW_SCALE_FACTORS'  : [1.4, 1.6, 1.8, 2.0],
                'VDW_POINT_DENSITY'  : int(gridPsi4),
                'RESP_A'             : 0.0005,
                'RESP_B'             : 0.1,
                'BASIS_ESP' : basisPsi4,
                'METHOD_ESP' : methodPsi4,
    }

    psi4.set_options({'reference': 'uhf'})
    charge = resp.resp([mol], options)[1]
    Logs.message("Done computing psi4 resp")
    return charge

def main():
    integrations = [enums.Integrations.hydrogen]
    plugin = nanome.Plugin('Hydrogens', 'A Nanome Plugin to add/remove hydrogens to/from structures', 'Hydrogens', True, integrations=integrations)
    plugin.set_plugin_class(Hydrogens)
    plugin.run()

if __name__ == '__main__':
    main()
