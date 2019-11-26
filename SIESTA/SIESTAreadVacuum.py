"""
This file reads the the file system_label.out and extracts/calculates data/values from it
"""

class SiestaReadOut():
    def __init__(self, system_label):
        self.system_label = system_label
        self.file = []
        f = open(f"{self.system_label}.out", "r")
        for line in f:
            self.file.append(line)
        f.close()

    def read_vacuum(self):
        """
        |This function sets the maximum and the mean values of the vacuum and the units. 
        |Inputs:
        |    None
        |Outputs:
        |    None
        """
        self.Vac_max = 0
        self.Vac_mean = 0
        self.Vac_units = 0
        for line in self.file:
            if "dhscf: Vacuum level (max, mean)" in line:
                x = line.strip("dhscf: Vacuum level (max, mean)=")
                x = x.split()
                self.Vac_max = float(x[0])
                self.Vac_mean = float(x[1])
                self.Vac_units = x[2]
                break

    def read_fermi(self):
        """
        |This function reads the fermi level. For now it is set only to read in eV
        |Inputs:
        |    None
        |Outputs:
        |    None
        """
        self.Ef = 0
        for line in self.file:
            if "siesta:         Fermi =" in line:
                x = line.strip("siesta:         Fermi =")
                x = x.split()
                self.Ef = float(x[0])
                break

    def calculate_work_function(self):
        """
        |This function calculates the work functions and saves in self.WF 
        |Inputs:
        |    None
        |Outputs:
        |    None
        """
        try:
            # If fermi and vacuum levels are calculated we can put just use them
            self.WF = self.Vac_mean - self.Ef
        except:
            # If fermi and vacuum levels are not calculated we need to calculate them first
            self.read_fermi()
            self.read_vacuum()
            self.WF = self.Vac_mean - self.Ef

    def read_total_energy(self):
        """
        |This function reads the total_energy. For now it is set only to read in eV
        |Inputs:
        |    None
        |Outputs:
        |    None
        """
        self.Ef = 0
        for line in self.file:
            if "siesta:         Total =" in line:
                x = line.strip("siesta:         Total =")
                x = x.split()
                self.E_total = float(x[0])
                break


if __name__ == "__main__":
    vac = SiestaReadOut("Fe")
    # vac.read_vacuum()
    # vac.read_fermi()
    vac.read_total_energy()
    print(vac.E_total)

    vac.calculate_work_function()
    print(vac.Ef)
    print(vac.Vac_mean)
    print(vac.WF)
