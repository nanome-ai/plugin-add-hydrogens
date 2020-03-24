import nanome
import tempfile
import shutil
import subprocess
from nanome.util import Logs

class Hydrogens(nanome.PluginInstance):
    def start(self):
        self.request = None
        self.temp_dir = tempfile.TemporaryDirectory()
        self.processes = []
        self.complex_input_files  = []
        self.serials = []
        self.matrices = []
        self.complex_output_files = []
        self.integration.hydrogen_add = self.add_H

    def update(self):
        if self.request and self.check_processes():
            complexes = []
            for i, complex_file in enumerate(self.complex_output_files):
                complex = nanome.structure.Complex.io.from_pdb(path=complex_file.name)
                complex.index = self.serials[i][0]
                for j, molecule in enumerate(complex.molecules):
                    molecule.index = self.serials[i][j]
                self.convert_to_relative_position(complex, i)
                complexes.append(complex)
                for atom in list(complex.atoms):
                    Logs.debug(f'{atom.symbol}: {atom.position}')
            self.request.send_response(complexes)
            self.clear()
        
    def clear(self):
        self.request = None
        self.complex_input_files.clear()
        self.complex_output_files.clear()
        self.serials.clear()
        self.matrices.clear()

    def add_H(self, request):
        Logs.debug('doing things!')
        self.request = request
        complexes = request.get_args()
        self.add_Hs_nanobabel(complexes)

    def on_stop(self):
        shutil.rmtree(self.temp_dir.name)

    def add_Hs_nanobabel(self, complexes):
        for complex in complexes:
            self.add_matrices_for_complex(complex)
            for atom in list(complex.atoms):
                Logs.debug(f'{atom.symbol}: {atom.position}')
            self.serials.append([complex.index])
            for molecule in complex.molecules:
                self.serials[-1].append(molecule.index)

            infile = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=self.temp_dir.name)
            outfile = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb", dir=self.temp_dir.name)
            self.complex_input_files.append(infile)
            self.complex_output_files.append(outfile)
            complex.io.to_pdb(infile.name)
            args = ['nanobabel', 'hydrogen', '-add', '-i', infile.name, '-o', outfile.name]
            self.processes.append(subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE))

    def check_processes(self):
        if len(self.processes) == 0:
            return True
        for process in self.processes:
            if process.poll() != None:
                index = self.processes.index(process)
                self.processes.remove(process)
        return False

    def add_matrices_for_complex(self, complex):
        self.matrices.append(complex.get_complex_to_workspace_matrix())

    def convert_to_relative_position(self, complex, i):
        complex.position = complex.get_workspace_to_complex_matrix() * complex.position
        complex.position = self.matrices[i] * complex.position

def main():
    plugin = nanome.Plugin('Hydrogens', 'A nanome integration plugin to add and remove hydrogens to/from structures', 'Hydrogens', False)
    plugin.set_plugin_class(Hydrogens)
    plugin.run('127.0.0.1', 8888)

if __name__ == '__main__':
    main()