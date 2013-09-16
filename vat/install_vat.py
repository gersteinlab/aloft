#!/usr/bin/env python
import subprocess
import os
import shutil

#this procedure extracts the .tar.gz files, creates a directory to install to, cd's into the source directory, configures it, builds it, and installs it
def installModule(moduleName, compileAsStatic=True, dependencies=[], sourceFilesToKill=[]):
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
	command = ["./configure", "--prefix=" + installPath] + staticFlags + cFlags
	os.system(" ".join(command)) #this command does not work very well when using subprocess.call, not sure why

	for sourceFileToKill in sourceFilesToKill:
		filepath = os.path.join("src", sourceFileToKill)
		if os.path.exists(filepath):
			os.remove(filepath)
		
		newFile = open(filepath, "w")

		newContents = "#include <stdio.h>\nint main(void) {\nprintf(\"%s is not compiled properly, get it from somewhere else\\n\");\nreturn 0;\n}" % (sourceFileToKill)
		newFile.write(newContents)
		newFile.close()

	subprocess.call(["make"])
	subprocess.call(["make", "install"])
	os.chdir(currentDirectory)
	shutil.rmtree(moduleName)

def buildVat():
	installModule('gsl-1.15')
	installModule('libgd-2.1.0')
	installModule('libbios-1.0.0', True, ['gsl-1.15'])
	installModule('libpng-1.6.4')
	installModule('freetype-2.4.0')
	#don't compile vcf2images.c as it causes potential trouble compiling
	#we'll first run configure then replace the source file with a dummy
	#this is a little hacky but automake doesn't make this very easy to do (i.e, can't just alter source file before running configure)
	installModule('vat-2.0.1', False, ['gsl-1.15', 'libbios-1.0.0', 'libgd-2.1.0', 'libpng-1.6.4', 'freetype-2.4.0'], ['vcf2images.c'])

if __name__ == "__main__":
	buildVat()
