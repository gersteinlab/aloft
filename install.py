#!/usr/bin/env python
import subprocess
import os
from vat.install_vat import buildVat

subprocess.call(['cc', 'gerprate.c', '-Wall', '-o', 'gerprate'])

os.chdir('vat')
buildVat()