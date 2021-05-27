import nanome
import tempfile
import shutil
import subprocess
from nanome.util import Logs
from nanome.util.enums import Integrations

PDBOptions = nanome.util.complex_save_options.PDBSaveOptions()
PDBOptions.write_bonds = True


class Hydrogens(nanome.PluginInstance):
    def start(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.processes = []
        self.integration.hydrogen_add = self.add_H
        self.integration.hydrogen_remove = self.rem_H
        self.complexes = None

    def on_run(self):
        def hydrogens(complexes):
            self.complexes = complexes
            req_rem = MemoryObject(
                {}, {'get_args': lambda: complexes, 'send_response': lambda complexes: self.add_H(req_add, complexes)})

            req_add = MemoryObject(
                {}, {'get_args': lambda: complexes, 'send_response': lambda complexes: self.update_structures_deep(complexes)})

            self.rem_H(req_rem)

        def get_selected(complexes):
            selected_complex_ids = [complex.index for complex in complexes if complex.get_selected()]

            if len(selected_complex_ids) == 0:
                self.plugin.send_notification(nanome.util.enums.NotificationTypes.warning, "Please select a complex.")
                return

            self.request_complexes(selected_complex_ids, hydrogens)

        self.request_complex_list(get_selected)

    def update(self):
        for process in reversed(self.processes):
            process[0].communicate()
            if process[0].poll() == None:
                continue

            complexes = []
            complex = nanome.structure.Complex.io.from_sdf(path=process[4].name)
            complex.position = self.complexes[0].position
            complex.rotation = self.complexes[0].rotation
            complex.index = process[1]
            for i, molecule in enumerate(complex.molecules):
                molecule.index = process[2][i]
            complexes.append(complex)
            process[5].send_response(complexes)
            Logs.debug('Removing process')
            self.processes.remove(process)

    def add_H(self, request, passed_complexes=None):
        Logs.debug('Add H')
        self.complexes = passed_complexes if passed_complexes is not None else request.get_args()
        self.exec_nanobabel('-add', self.complexes, request)

    def rem_H(self, request):
        Logs.debug('Remove H')
        self.complexes = request.get_args()
        self.exec_nanobabel('-del', self.complexes, request)

    def on_stop(self):
        shutil.rmtree(self.temp_dir.name)

    def exec_nanobabel(self, arg, complexes, request):
        for complex in complexes:
            mol_serials = []
            for molecule in complex.molecules:
                mol_serials.append(molecule.index)

            infile = tempfile.NamedTemporaryFile(
                delete=False, suffix=".pdb", dir=self.temp_dir.name)
            outfile = tempfile.NamedTemporaryFile(
                delete=False, suffix=".sdf", dir=self.temp_dir.name)
            complex.io.to_pdb(infile.name, PDBOptions)
            args = ['nanobabel', 'hydrogen', arg,
                    '-i', infile.name, '-o', outfile.name]
            p = subprocess.Popen(
                args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.processes.append(
                (p, complex.index, mol_serials, infile, outfile, request))

    def check_processes(self):
        if len(self.processes) == 0:
            return True
        for process in self.processes:
            if process[0].poll() != None:
                self.processes.remove(process)
        return False


def main():
    plugin = nanome.Plugin('Hydrogens', 'A nanome integration plugin to add and remove hydrogens to/from structures', 'Hydrogens', False, integrations=[Integrations.hydrogen])
    plugin.set_plugin_class(Hydrogens)
    plugin.run()


if __name__ == '__main__':
    main()


class MemoryObject:
    def __init__(self, obj, functions):
        self.__dict__.update(obj)
        for function in functions:
            self.__setattr__(function, functions[function])
