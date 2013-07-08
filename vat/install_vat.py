#!/usr/bin/env python
import subprocess
import os
import shutil

#this procedure extracts the .tar.gz files, creates a directory to install to, cd's into the source directory, configures it, builds it, and installs it
def installModule(moduleName, dependencies=[]):
	if os.path.exists(moduleName):
		shutil.rmtree(moduleName)
	
	subprocess.call(['tar', '-zxvf', moduleName+".tar.gz"])
	currentDirectory = os.getcwd()
	installDirectory = moduleName+"-install"
	
	if os.path.exists(installDirectory):
		shutil.rmtree(installDirectory)
	
	os.mkdir(installDirectory)

	installPath = os.path.join(currentDirectory, installDirectory)
	os.chdir(moduleName)

	dependencePrefixes = ["-I%s -L%s" % (os.path.join(os.path.join(currentDirectory, dependence+"-install"), "include"), os.path.join(os.path.join(currentDirectory, dependence+"-install"), "lib")) for dependence in dependencies]
	if len(dependencePrefixes) > 0:
		cFlags = ['CFLAGS="%s"' % (" ".join(dependencePrefixes))]
	else:
		cFlags = []
	command = ["./configure", "--prefix=" + installPath, "--enable-static", "--disable-shared"] + cFlags
	os.system(" ".join(command)) #this command does not work very well when using subprocess.call, not sure why
	subprocess.call(["make"])
	subprocess.call(["make", "install"])
	os.chdir(currentDirectory)
	shutil.rmtree(moduleName)

def buildVat():
	installModule('gsl-1.15')
	installModule('libgd-2.1.0')
	installModule('libbios-1.0.0', ['gsl-1.15'])
	installModule('vat-2.0.1', ['gsl-1.15', 'libbios-1.0.0', 'libgd-2.1.0'])

if __name__ == "__main__":
	buildVat()
