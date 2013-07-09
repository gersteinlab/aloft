#!/usr/bin/env python
import subprocess, os, shutil
from vat.install_vat import buildVat

subprocess.call(['cc', 'gerprate.c', '-Wall', '-o', 'gerprate'])

networkModule = 'networkx'
networkDirectory = networkModule+'-1.7'

if os.path.exists(networkDirectory):
	shutil.rmtree(networkDirectory)
subprocess.call(['tar', '-zxvf', networkDirectory+".tar.gz"])

if os.path.exists(networkModule):
	shutil.rmtree(networkModule)
shutil.move(os.path.join(networkDirectory, networkModule), networkModule)

shutil.rmtree(networkDirectory)

os.chdir('vat')
buildVat()