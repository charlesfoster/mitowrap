#!/usr/bin/env python
import sys
import os 

with open(sys.argv[1],'r') as f:
    using_conda = f.read().strip()

if using_conda == "True":
    print("Installing ete3 NCBI data")
    from ete3 import NCBITaxa
    ncbi = NCBITaxa()
    os.system("echo 'data added' > {}".format(sys.argv[2]))
else:
    os.system("echo 'data not needed to be added' >  {}".format(sys.argv[2]))
