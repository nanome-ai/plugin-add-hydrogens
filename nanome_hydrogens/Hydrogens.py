import nanome
import tempfile
from nanome.util import async_callback, Logs, Process
from nanome.util.enums import Integrations, NotificationTypes
from .ComplexUtils import ComplexUtils as utils

PDBOptions = nanome.util.complex_save_options.PDBSaveOptions()
PDBOptions.write_bonds = True

class Hydrogens(nanome.AsyncPluginInstance):
    def start(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.input_file = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb', dir=self.temp_dir.name)
        self.output_file = tempfile.NamedTemporaryFile(delete=False, suffix='.sdf', dir=self.temp_dir.name)
        self.integration.hydrogen_add = self.add_hydrogens
        self.integration.hydrogen_remove = self.remove_hydrogens

    @async_callback
    async def on_run(self):
        # get selected
        shallow = await self.request_complex_list()
        indices_selected = [c.index for c in shallow if c.get_selected()]

        # were there any?
        if not indices_selected:
            self.send_notification(NotificationTypes.warning, 'Please select a complex.')
            return

        # remove hydrogens
        deep = await self.request_complexes(indices_selected)
        deep = await self.remove_hydrogens(complexes=deep)

        # stop if something's wrong
        if not deep:
            return

        # add hydrogens
        await self.add_hydrogens(complexes=deep, upload=True)

    def on_stop(self):
        self.temp_dir.cleanup()

    @async_callback
    async def add_hydrogens(self, request=None, complexes=None, upload=False):
        Logs.debug('Add H')
        new_complexes = await self.run_hydrogens('-add', request, complexes, upload, 'Could not add hydrogens.')
        return new_complexes

    @async_callback
    async def remove_hydrogens(self, request=None, complexes=None, upload=False):
        Logs.debug('Remove H')
        new_complexes = await self.run_hydrogens('-del', request, complexes, upload, 'Could not remove hydrogens.')
        return new_complexes

    @async_callback
    async def run_hydrogens(self, arg, request, complexes, upload=False, err_message='An error occurred'):
        # Run or Integration?
        to_update = request.get_args() if request else complexes

        for ci, complex in enumerate(to_update):
            complex.io.to_pdb(self.input_file.name, PDBOptions)

            p = Process()
            p.executable_path = 'nanobabel'
            p.args = ['hydrogen', arg, '-i', self.input_file.name, '-o', self.output_file.name]
            p.output_text = True
            p.on_error = Logs.error
            p.on_output = Logs.debug

            # if error, notify user and return
            exit_code = await p.start()
            if exit_code:
                self.send_notification(NotificationTypes.error, err_message)
                Logs.error(err_message)
                return

            updated_complex = nanome.api.structure.Complex.io.from_sdf(path=self.output_file.name)
            molecules = list(complex.molecules)
            for mi, updated_molecule in enumerate(updated_complex.molecules):
                updated_molecule.index = molecules[mi].index
            # if polar setting:
            utils.markPolarHydrogens(updated_complex)
            utils.reidentify(updated_complex, complex)
            to_update[ci] = updated_complex

        if request:
            request.send_response(to_update)
        elif upload:
            self.update_structures_deep(to_update)

        return to_update

def main():
    integrations = [Integrations.hydrogen]
    plugin = nanome.Plugin('Hydrogens', 'A Nanome Plugin to add/remove hydrogens to/from structures', 'Hydrogens', False, integrations=integrations)
    plugin.set_plugin_class(Hydrogens)
    plugin.run()

if __name__ == '__main__':
    main()
