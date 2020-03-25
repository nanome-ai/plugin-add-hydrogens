import nanome
import tempfile
import shutil
import subprocess
from nanome.util import Logs

class Hydrogens(nanome.PluginInstance):
    def start(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.processes = []
        self.integration.hydrogen_add = self.add_H
        self.integration.hydrogen_remove = self.rem_H

    def update(self):
        for process in reversed(self.processes):
            process[0].communicate()
            if process[0].poll() == None:
                continue

            complexes = []
            complex = nanome.structure.Complex.io.from_sdf(path=process[4].name)
            complex.index = process[1]
            for i, molecule in enumerate(complex.molecules):
                molecule.index = process[2][i]
            complexes.append(complex)
            process[5].send_response(complexes)
            Logs.debug('Removing process')
            self.processes.remove(process)

    def add_H(self, request):
        Logs.debug('Add H')
        complexes = request.get_args()
        self.exec_nanobabel('-add', complexes, request)

    def rem_H(self, request):
        Logs.debug('Remove H')
        complexes = request.get_args()
        self.exec_nanobabel('-del', complexes, request)

    def on_stop(self):
        shutil.rmtree(self.temp_dir.name)

    def exec_nanobabel(self, arg, complexes, request):
        for complex in complexes:
            mol_serials = []
            for molecule in complex.molecules:
                mol_serials.append(molecule.index)

            infile = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=self.temp_dir.name)
            outfile = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf", dir=self.temp_dir.name)
            complex.io.to_pdb(infile.name)
            args = ['nanobabel', 'hydrogen', arg, '-i', infile.name, '-o', outfile.name]
            p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.processes.append((p, complex.index, mol_serials, infile, outfile, request))

    def check_processes(self):
        if len(self.processes) == 0:
            return True
        for process in self.processes:
            if process[0].poll() != None:
                self.processes.remove(process)
        return False

def main():
    plugin = nanome.Plugin('Hydrogens', 'A nanome integration plugin to add and remove hydrogens to/from structures', 'Hydrogens', False)
    plugin.set_plugin_class(Hydrogens)
    plugin.run('127.0.0.1', 8888)

if __name__ == '__main__':
    main()