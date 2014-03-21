#!/usr/bin/env python
import subprocess
import os
import shutil

GSL_PATH = 'gsl-1.15'
BIOS_PATH = 'libbios-1.0.0'
JPEG_PATH = 'jpeg-9'
PNG_PATH = 'libpng-1.6.4'
GD_PATH = 'libgd-2.1.0'
FREETYPE_PATH = 'freetype-2.4.0'

#this procedure extracts the .tar.gz files, creates a directory to install to, cd's into the source directory, configures it, builds it, and installs it
def installModule(moduleName, compileAsStatic=True, dependencies=[], extraFlags=[]):
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

	if compileAsStatic:
		staticFlags = ["--enable-static", "--disable-shared"]
	else:
		staticFlags = []

	command = ["./configure", "--prefix=" + installPath] + staticFlags + extraFlags + cFlags
	os.system(" ".join(command)) #this command does not work very well when using subprocess.call, not sure why

	subprocess.call(["make"])
	subprocess.call(["make", "install"])
	os.chdir(currentDirectory)
	shutil.rmtree(moduleName)

	return installPath

def buildVat():
	installModule(GSL_PATH)
	installModule(BIOS_PATH, True, [GSL_PATH])

	jpegInstallPath = installModule(JPEG_PATH)
	pngInstallPath = installModule(PNG_PATH)

	installModule(GD_PATH, True, [], ["--with-jpeg=%s" % jpegInstallPath, "--with-png=%s" % pngInstallPath])
	installModule(FREETYPE_PATH)

	vatInstallPath = installModule('vat-2.0.1', False, [GSL_PATH, JPEG_PATH, PNG_PATH, BIOS_PATH, GD_PATH, FREETYPE_PATH])

	return os.path.join(vatInstallPath, 'bin')

if __name__ == "__main__":
	buildVat()
