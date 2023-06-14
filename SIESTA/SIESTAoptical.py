import matplotlib.pyplot as plt

class optical:
    def __init__(self) -> None:
        self.file_name = "Optical"
        self.set_x_range = True
        self.xlim_low = 0
        self.xlim_high = 10
        self.kwargs = []  # for all the matplotlib pltotting kwargs
        self.plt_show = True
        self.labels = []
        self.e1 = []
        self.e2 = []
        self.abs = []
        self.fig_extension = "pdf"

    def load(self, out_file, label=None, **kwargs) -> None:
        self.kwargs.append(kwargs)
        if label == None: self.labels.append(out_file.SystemLabel)
        else: self.labels.append(label)
        file = open(f"{out_file.folder_path}/epsilon_real.out", 'r')
        data = [line for line in file]
        file.close()
        energies = []
        vals = []
        for line in data:
            energies.append(float(line.split()[0]))
            vals.append(float(line.split()[1]))
        self.e1.append([energies, vals])

        file = open(f"{out_file.folder_path}/epsilon_img.out", 'r')
        data = [line for line in file]
        file.close()
        energies = []
        vals = []
        for line in data:
            energies.append(float(line.split()[0]))
            vals.append(float(line.split()[1]))
        self.e2.append([energies, vals])
        
        file = open(f"{out_file.folder_path}/absorp_coef.out", 'r')
        data = [line for line in file]
        file.close()
        energies = []
        vals = []
        for line in data:
            energies.append(float(line.split()[0]))
            vals.append(float(line.split()[1]))
        self.abs.append([energies, vals])

    def plot_e2(self) -> None:
        for i,val in enumerate(self.e2):
            plt.plot(val[0], val[1], label = f"{self.labels[i]}", **self.kwargs[i])
        plt.xlabel(f'Energy (eV)')
        plt.ylabel(f'Im[$\epsilon$]')
        plt.legend(loc='upper right')
        # plt.title(f"{self.plt_title}")
        if self.set_x_range == True:
            plt.xlim([self.xlim_low,self.xlim_high])
        # if self.set_y_range == True:
        #     plt.ylim([self.ylim_low,self.ylim_high])
        plt.tight_layout()
        plt.savefig(f"{self.file_name}_e2.{self.fig_extension}")
        if self.plt_show: plt.show()
        plt.close()

    def plot_absorption(self, ylog = False) -> None:
        for i,val in enumerate(self.abs):
            plt.plot(val[0], val[1], label = f"{self.labels[i]}", **self.kwargs[i])
        plt.xlabel(f'Energy (eV)')
        plt.ylabel(f'Absorption (cm$^{{-1}}$)')
        if ylog: plt.yscale("log")
        plt.legend(loc='upper right')
        # plt.title(f"{self.plt_title}")
        if self.set_x_range == True:
            plt.xlim([self.xlim_low,self.xlim_high])
            plt.ylabel(f'Absorption (log) (cm$^{{-1}}$)')
        # if self.set_y_range == True:
        #     plt.ylim([self.ylim_low,self.ylim_high])
        plt.tight_layout()
        plt.savefig(f"{self.file_name}_absorption.{self.fig_extension}")
        if self.plt_show: plt.show()
        plt.close()
        