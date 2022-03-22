import os
import commands
import sys

#os.system('source /opt/app-sync/intel/bin/iccvars.sh intel64')

try:
 maxrun = int(sys.argv[1])
except:
 maxrun = 3

try:
    machine = sys.argv[2]
except:
    machine = "meson"

if machine == "meson" or machine == "Meson" or machine == "MESON":
    s = commands.getstatusoutput("ps uaxw |grep OBCZRH |grep -v -c grep")
    running = int(s[1])
else:
    running = 0

free = maxrun - running

if free > 0:
 files = [f for f in os.listdir(os.path.expanduser('~/Documents/OBCZRH/Jobs'))]

 n = 0
 while n < len(files) and n < free:
  #print(files[n])
  thisfile = os.path.expanduser('~/Documents/OBCZRH/Jobs/./' + files[n])
  os.system(thisfile)
  n += 1
