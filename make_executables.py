import shutil, os

scripts = ['aloft.py', 'vcf_sort.py']
for script in scripts:
	executableName = script.rstrip(".py")
	shutil.copy2(script, executableName)
	os.chmod(executableName, os.stat(executableName).st_mode | 0o111) #set executable bits
