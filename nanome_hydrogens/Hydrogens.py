import nanome
import tempfile
import shutil
from nanome.util import Logs

class Hydrogens(nanome.PluginInstance):
    def start(self):
        self.request = None
        self.temp_dir = tempfile.TemporaryDirectory()
        self.processes = []
        self.complex_input_files  = []
        self.complex_serials = []
        self.molecule_serials = []
        self.complex_output_files = []
        self.integration.hydrogen_add = self.add_H

    def update(self):
        if self.request and self.check_processes():
            complexes = []
            for i, complex_file in enumerate(self.complex_output_files):
                complex = nanome.structure.Complex.from_pdb(complex_file.name)
                complex.index = self.complex_serials[i]
                for i, molecule in enumerate(complex.get_molecules()):
                    molecule.index = self.molecule_serials[i]
                complexes.append(complex)
            self.request.send_response(complexes)
            self.request = None
            self.complex_serials = []
            self.molecule_serials = []

    def add_H(self, request):
        Logs.debug('doing things!')
        self.request = request
        complexes = request.get_args()
        self.add_Hs_nanobabel(complexes)

    def on_stop(self):
        shutil.rmtree(self.temp_dir.name)

    def add_Hs_nanobabel(self, complexes):
        if len(self._complex_files):
            return

        for complex in complexes:
            self.complex_serials.append(complex.index)
            for molecule in complex.get_molecules():
                self.molecule_serials.append(molecule.index)

            infile = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=self.temp_dir.name)
            outfile = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=self.temp_dir.name)
            self.complex_input_files.append(infile)
            self.complex_output_files.append(outfile)
            complex.io.to_pdb(infile.name)
            args = ['nanobabel', 'hydrogen', '-add', '-i', infile.name, '-o', self.outfile.name]
            self.processes.append(subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE))

    def check_processes(self):
        if len(self.processes) == 0:
            return True
        for process in self.processes:
            if self._smina_process.poll() != None:
                self.processes.remove(process)
        return False

def main():
    plugin = nanome.Plugin('Hydrogens', 'A nanome integration plugin to add and remove hydrogens to/from structures', 'Hydrogens', False)
    plugin.set_plugin_class(Hydrogens)
    plugin.run('127.0.0.1', 8888)

if __name__ == '__main__':
    main()