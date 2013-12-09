import shutil, os

shutil.copy2('aloft.py', 'aloft')
os.chmod("aloft", os.stat('aloft').st_mode | 0o111) #set executable bits
