#!/usr/bin/env python

import os, subprocess, shutil
from vat.install import buildVat

def getDirectory():
	return os.path.dirname(os.path.realpath(__file__))

def installVat():
	currentDirectory = getDirectory()
	os.chdir(os.path.join(currentDirectory, 'vat'))
	vatBinPath = buildVat()
	os.chdir(currentDirectory)

	newBinPath = 'vat-bin'
	if os.path.exists(newBinPath):
		shutil.rmtree(newBinPath)
		
	subprocess.call(['cp', '-r', vatBinPath, newBinPath])

if __name__ == "__main__":
	installVat()