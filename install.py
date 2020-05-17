"""This file will create contains installation scripts"""
import logging

logging.basicConfig(filename='installation_log.txt', level=logging.DEBUG, format='\
%(asctime)s - %(levelname)s - %(message)s')

class Install():
    def __init__(self, source_folder):
        self.source_folder = source_folder

    def check(self):
        pass

    def install(self):
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
        pass

    def check(self):
        pass

    def install(self):
        pass

if __name__ == "__main__":
    pass