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

class Install():
    def __init__(self):
        self.source_folder = os.getcwd()
        self.system_version = get_system_version()

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
        return answer

    def check_master(self):
        """
        This is the master checking method specific other checks should be done with the other functions in the daughter class methods
        """
        check_var = []
        check_var.append(self.check_system_version())
        if False in chark_var: return False
        else: return True

    def check_system_version(self):
        """
        This method returns a boolean True if all checks are postive and therefore
        intallatin can proceed. Else will return false
        """
        if self.system_version == "win32":
            is_wls = subprocess.Popen("wsl", shell=True, stdout=subprocess.PIPE)
            answer = is_wls.stdout.read()
            print(answer)
            logging.critical("This application needs a unix like system to work")
            logging.critical("Seems like wsl is not installed - Quitting installation")
            return False
        else:
            return True

    def install(self):
        """
        This method is the main install method
        """
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

    def check_this(self):
        pass

    def _install(self):
        pass

class InstallQuantumEspresso(Install):
    def __init__(self):
        self.installation_package = "Quantum Espresso"
        super().__init__(self)

    def check_this(self):
        pass

    def _install(self):
        pass

if __name__ == "__main__":
    pass