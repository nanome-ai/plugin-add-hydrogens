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
        self.input_file = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=self.temp_dir.name)
        self.output_file = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf", dir=self.temp_dir.name)
        self.integration.hydrogen_add = self.add_hydrogens
        self.integration.hydrogen_remove = self.remove_hydrogens

    @async_callback
    async def on_run(self):
        # get selected
        shallow = await self.request_complex_list()
        indices_selected = [c.index for c in shallow if c.get_selected()]

        # were there any?
        if not indices_selected:
            self.send_notification(NotificationTypes.warning, "Please select a complex.")
            return

        # remove hydrogens
        deep = await self.request_complexes(indices_selected)
        deep = await self.remove_hydrogens(None, deep, upload=False)

        # stop if something's wrong
        if not deep:
            return

        # add hydrogens
        await self.add_hydrogens(None, deep, upload=True)

    def on_stop(self):
        self.temp_dir.cleanup()

    @async_callback
    async def add_hydrogens(self, request, complexes=None, upload=False):
        Logs.debug('Add H')
        new_complexes = await self.run_hydrogens('-add', request, complexes, upload, err_message="Could not add hydrogens.")
        return new_complexes

    @async_callback
    async def remove_hydrogens(self, request, complexes=None, upload=False):
        Logs.debug('Remove H')
        new_complexes = await self.run_hydrogens('-del', request, complexes, upload, err_message="Could not delete hydrogens.")
        return new_complexes

    @async_callback
    async def run_hydrogens(self, arg, request, complexes, upload=False, err_message='An error occured'):
        # Run or Integration?
        to_update = request.get_args() if request else complexes

        for ci, complex in enumerate(to_update):
            
            # nanobabel as Process API Process
            complex.io.to_pdb(self.input_file.name, PDBOptions)
            p = create_process('nanobabel', ['hydrogen', arg, '-i', self.input_file.name, '-o', self.output_file.name])
            
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

def create_process(path, args):
    p = Process()
    p.executable_path = path
    p.args = args
    p.output_text = True
    p.on_error = Logs.error
    p.on_output = Logs.warning
    p.on_done = Logs.message
    return p

def main():
    plugin = nanome.Plugin('(Polar) Hydrogens', 'A nanome integration plugin to add and remove hydrogens to/from structures', 'Hydrogens', False, integrations=[Integrations.hydrogen])
    plugin.set_plugin_class(Hydrogens)
    plugin.run()


if __name__ == '__main__':
    main()


class MemoryObject:
    def __init__(self, obj, functions):
        self.__dict__.update(obj)
        for function in functions:
            self.__setattr__(function, functions[function])
