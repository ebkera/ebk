"""This file will create contains installation scripts"""
import sys
import logging

# Some global settings here
logging.basicConfig(filename='installation.log.era', level=logging.DEBUG, format='\
%(asctime)s - %(levelname)s - %(message)s')

def get_system_version():
    version = sys.platform
    logging.info(f"System platform detected: {version}")
    return version

class Install():
    def __init__(self, source_folder):
        self.source_folder = source_folder
        self.system_version = get_system_version()

    def check(self):
        pass

    def install(self):
        pass

class win32(Install):
    def __init__(self):
        pass

class InstallSIESTA(Install):
    def __init__(self, source_folder):
        super().__init__(source_folder)

    def check(self):
        pass

    def install(self):
        pass

class InstallSIESTA(Install):
    def __init__(self):
        super().__init__(source_folder)

    def check(self):
        pass

    def install(self):
        pass

if __name__ == "__main__":
    pass