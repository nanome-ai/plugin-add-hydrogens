import nanome
from nanome.util import async_callback, Logs
from nanome.util.enums import Integrations, NotificationTypes
from rdkit import Chem
import shutil
import tempfile

PDBOptions = nanome.util.complex_save_options.PDBSaveOptions()
PDBOptions.write_bonds = True

class Hydrogens(nanome.PluginInstance):
    def start(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.input_file = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=self.temp_dir.name)
        self.output_file = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=self.temp_dir.name)
        self.processes = []
        self.integration.hydrogen_add = self.add_H
        self.integration.hydrogen_remove = self.rem_H

    @async_callback
    async def on_run(self):
        shallow = await self.request_complex_list(get_selected)
        # get selected complexes
        selected_complex_ids = [complex.index for complex in shallow if complex.get_selected()]
        if not selected_complex_ids:
            self.send_notification(NotificationTypes.warning, "Please select a complex.")
            return

        deep = await self.request_complexes(selected_complex_ids)
        self.add_H(self, complexes=deep, upload=False)
        self.rem_H(self, complexes=deep, upload=True)


    def add_H(self, request, complexes=None, upload=False):
        Logs.debug('Add Hs')
        to_be_Hd = complexes if complexes else request.get_args()
        for complex in to_be_Hd:
            # nanome complex -> pdb -> rdkit molecule
            pdb = complex.io.to_pdb(self.input_file, PDBOptions)
            rdmol = Chem.rdmolfiles.MolFromPDBFile(pdb)
            # add Hs to rdkit molecule
            Chem.AddHs(rdmol)
            Chem.rdmolfiles.MolToPDBFile(mol, self.output_file.name)
            # rdkit molecule -> pdb -> nanome complex
            H_complex = nanome.Complex.io.from_pdb(self.output_file.name)
            H_complex.index = complex.index

            # upload
            if upload:
                self.update_structures_deep([H_complex])

    def rem_H(self, request, complexes=None, upload=False):
        Logs.debug('Remove Hs')
        to_be_unHd = complexes if complexes else request.get_args()
        for complex in self.to_be_unHd:
            # nanome complex -> pdb -(remove Hs)-> rdkit molecule
            pdb = complex.io.to_pdb(self.input_file, PDBOptions)
            rdmol = Chem.rdmolfiles.MolFromPDBFile(pdb, removeHs=True)
            # rdkit molecule -> pdb -> nanome complex
            Chem.rdmolfiles.MolToPDBFile(rdmol, self.output_file.name)
            NH_complex = nanome.Complex.io.from_pdb(self.output_file.name)
            NH_complex.index = complex.index

            # upload
            if upload:
                self.update_structures_deep([NH_complex])

    def on_stop(self):
        shutil.rmtree(self.temp_dir.name)

def main():
    plugin = nanome.Plugin('Hydrogens', 'A nanome integration plugin to add and remove hydrogens to/from structures', 'Hydrogens', False, integrations=[Integrations.hydrogen])
    plugin.set_plugin_class(Hydrogens)
    plugin.run()


if __name__ == '__main__':
    main()