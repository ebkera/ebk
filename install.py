"""This file will create contains installation scripts"""
import sys
import logging
import subprocess
import os

# Some global settings here
logging.basicConfig(filename='installation.log.era', level=logging.DEBUG, format='\
%(asctime)s - %(levelname)s - %(message)s')

def get_system_version():
    """Gets the system version."""
    return sys.platform

def install_openmpi():
    """
    Attempts to install Open MPI
    """
    pass

class Install():
    def __init__(self):
        self.source_folder = os.getcwd()
        self.system_version = get_system_version()
        self.parallel_installation = True

    def system_call(self, command):
        """
        This command calls the shell and does it system specifically. 
            If linux machine calls shell
            if windows then calls shell commands in wsl
        """
        if self.system_version == "win32":
            call = subprocess.Popen(f"wsl {command}", shell=True, stdout=subprocess.PIPE)
        else:
            call = subprocess.Popen(f"{command}", shell=True, stdout=subprocess.PIPE)
        answer = call.stdout.read()
        return str(answer)

    def update_packages(self):
        print("Updating system first")
        self.system_call("sudo apt-get update")
        self.system_call("sudo apt-get upgrade")

    def check_master(self):
        """
        This is the master checking method specific other checks should be done with the other functions in the daughter class methods.
        How the method works:
            It first calls the methods that has to be checked and gets boolean values that is saves into a list
            Then specific checks are performed by calling check_daughter method of the instance
            If they also return true then the installation can proceed
        """
        check_var = []
        # Updating system first
        self.update_packages()
        # First checking the master settings and installations
        check_var.append(self.check_system_version())
        if self.parallel_installation: check_var.append(self.check_open_mpi())
        # Then checking daughter specific settings and installations
        check_var.extend(self.check_daughter())
        if False in check_var: return False
        else: return True

    def check_system_version(self):
        """
        This method returns a boolean True if all checks are postive and therefore
        intallation can proceed. Else will return false
        """
        if self.system_version == "win32":
            is_wls = subprocess.Popen("wsl ls", shell=True, stdout=subprocess.PIPE)
            answer = str(is_wls.stdout.read())
            logging.warning("This application needs a unix like system to work")
            if "log" in answer:
                logging.info("Seems like wsl is installed - Continuing installation")
                return True
            else:
                logging.critical("Seems like wsl is not installed - Quitting installation (continuing with other checks)")
                print("Installation stopped: Check to see if wsl is installed (continuing with other checks)")
                return False
        else:
            return True

    def check_open_mpi(self):
        answer = self.system_call("mpirun --version")
        print(answer)
        logging.warning("Settings are set to use Parallel installation")
        if "MPI" in answer:
            logging.info("Seems like Open MPI is installed - Continuing installation")
            return True
        else:
            logging.warning("Seems like Open MPI is not installed - Attempting to install Open MPI")
            install_openmpi()
            answer = subprocess.Popen("wsl mpirun --version", shell=True, stdout=subprocess.PIPE)
            answer = str(answer.stdout.read())
            logging.info("Settings are set to use Parallel installation")
            if "MPI" in answer:
                logging.warning("Successfully installed Open MPI is installed - Continuing installation")
                return True
            else:
                logging.critical("Seems like Open MPI has not installed - Quitting installation (continuing with other checks)")
                print("Installation stopped: Check to see if OpenMPI is installed (continuing with other checks)")
                return False

    def install(self):
        """
        This method is the main install method
        """
        logging.info("***********************Installation Log************************")
        logging.info(f"System platform detected: {self.system_version}")
        logging.info(f"Starting installation of {self.installation_package}")
        test = self.check_master()
        print(test)


class win32(Install):
    def __init__(self):
        pass

class InstallSIESTA(Install):
    def __init__(self):
        self.installation_package = "SIESTA"
        super().__init__()

    def check_daughter(self):
        check_var = []


        return check_var

    def check_this(self):
        pass

    def _install(self):
        pass

class InstallQuantumEspresso(Install):
    def __init__(self):
        self.installation_package = "Quantum Espresso"
        super().__init__(self)

    def check_daughter(self):
        pass

    def check_this(self):
        pass

    def _install(self):
        pass

if __name__ == "__main__":
    pass