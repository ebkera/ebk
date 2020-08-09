from ase.io import read, write
from ase.atom import Atom
from ebk import get_machine_paths

xyz_path = get_machine_paths()["xyz"]
EDT12_path = f"{xyz_path}/1,2-ethaneDithiol_relaxed.xyz"
BDT12_path = f"{xyz_path}/1,2-benzeneDithiol_relaxed.xyz"
BDT14_path = f"{xyz_path}/1,4-benzeneDithiol_relaxed_KE50.xyz"

class Insert_ligand():
    """
    This class is made to be strictly for ligands so that we dont have to go about attaching ligands manually. 
    This class will isert a ligand at a given position and at a given orientation.
    """
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

    def attach(self, attach_to, attach_at, attach_through)
    """
    This is the main function that attacheds the ligand to a atoms type object
    attach_to: The atoms type object that the ligand will attach to
    attach_at: Will attach at this atom in the atoms object(attach_to) and if there is an atom there already it will remove it and attach the ligand
    attach_through: The ligand will remove this atom and the resulting will attach throught the resulting dangling bond. to the attach_at site in the ligand.
    """
        pass

class BDT14(Insert_ligand):
    def __init__(self, *args, **kwargs):
        self.name = kwargs.get("name", "BDT14")
        self.atoms = read(BDT14_path, index=None, format="xyz")
        super().__init__()

class BDT12(Insert_ligand):
    def __init__(self, *args, **kwargs):
        self.name = kwargs.get("name", "BDT12")
        self.atoms = read(BDT12_path, index=None, format="xyz")
        super().__init__()

class EDT12(Insert_ligand):
    def __init__(self, *args, **kwargs):
        self.name = kwargs.get("name", "EDT12")
        self.atoms = read(EDT12_path, index=None, format="xyz")
        super().__init__()