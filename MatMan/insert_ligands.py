from ase.io import read, write
from ase.atom import Atom
from ebk import get_machine_paths

xyz_path = get_machine_paths()["xyz"]
EDT_path = f"{xyz_path}/1,2-ethaneDithiol_relaxed.xyz"
BDT12_path = f"{xyz_path}/1,2-benzeneDithiol_relaxed.xyz"
BDT14_path = f"{xyz_path}/1,4-benzeneDithiol_relaxed_KE50.xyz"
print(EDT_path, BDT12_path, BDT14_path)

class Insert_ligand():
    def __init__(self, *args, **kwargs):
        pass

    def orient(self, direction):
        """
        Plan to make this into a function that takes in  direction and orients the ligands in that direction.
        This might have to be a ligands specific fucntion though.
        """
        if direction == [0,0,1]:
            self.atoms.rotate(90, 'y')
        if direction == [1,1,0]:
            self.atoms.rotate(45, 'z')

    def edit(self):
        self.atoms.edit()

    def get_chemical_symbols(self):
        return self.atoms.get_chemical_symbols()

    def get_positions(self):
        return self.atoms.get_positions()

    def list_all_atoms(self):
        syms = self.get_chemical_symbols()
        coor = self.get_positions()
        return [f"{syms[n]}: {coor[n]}" for n in range(len(syms))]

    def save(self, format):
        write(f"{self.name}.{format}", self.atoms)

class BDT14(Insert_ligand):
    def __init__(self, *args, **kwargs):
        self.name = kwargs.get("name", "BDT14")
        self.atoms = read(BDT14_path, index=None, format="xyz")
        super().__init__()

class BDT12(Insert_ligand):
    def __init__(self, *args, **kwargs):
        self.name = kwargs.get("name", "BDT12")
        self.atoms = read(BDT14_path, index=None, format="xyz")
        super().__init__()

class EDT12(Insert_ligand):
    def __init__(self, *args, **kwargs):
        self.name = kwargs.get("name", "EDT12s")
        self.atoms = read(BDT14_path, index=None, format="xyz")
        super().__init__()