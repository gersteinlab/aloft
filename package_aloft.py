import shutil, os, subprocess

#Make relevant scripts executable
scripts = ['aloft.py']
for script in scripts:
	executableName = script.rstrip(".py")
	shutil.copy2(script, executableName)
	os.chmod(executableName, os.stat(executableName).st_mode | 0o111) #set executable bits

#Set up networkx
networkModule = 'networkx'
networkDirectory = networkModule+'-1.8.1'

if os.path.exists(networkDirectory):
	shutil.rmtree(networkDirectory)

subprocess.call(['tar', '-zxvf', networkDirectory+".tar.gz"])

if os.path.exists(networkModule):
	shutil.rmtree(networkModule)

shutil.move(os.path.join(networkDirectory, networkModule), networkModule)

shutil.rmtree(networkDirectory)
