#!/usr/bin/env python
import sys
import os
#print sys.argv
#print len(sys.argv)
#print sys.argv[1]
if(len(sys.argv)> 1 and sys.argv[1]!= ''):
   print 'EXECUTING'
   fnm = sys.argv[1]
   print fnm
   #exit()
   for line in file(fnm):
	print line
   	os.system(line)


