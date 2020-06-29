from ase.io import read, write
from ase.atom import Atom

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

    def save(self, format):
        write(f"{self.name}.{format}", self.atoms)

class BDT14(Insert_ligand):
    def __init__(self, *args, **kwargs):
        self.name = kwargs.get("name", "BDT14")
        self.atoms = read("1,4-benzeneDithiol_relaxed_KE50.xyz", index=None, format="xyz")
        super().__init__()