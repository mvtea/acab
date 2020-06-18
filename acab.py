#-----------------------------------------------------------
# FILE: acab.py
# AUTHOR/UPDATED: Mason V. Tea / 6.16.20
# PURPOSE: Automated Characterization of Accreting Binaries.
#-----------------------------------------------------------

# Import libraries
from acab_funcs import fitting,ftest,log_info,choose_model
import numpy as np
from sherpa.astro.ui import *
import logging
import time

# Set start time
start_time = time.time()

# Input sources from file; If only one source, list of lists becomes list... source has 3 values associated
sources = np.loadtxt('srcs_2000_rk.csv', delimiter=',', unpack=False, dtype=str, comments='#')
if len(sources) == 3:
    sources = [sources]

# TEMPORARILY CONTRAIN NUMBER OF SOURCES FOR TESTING
# sources = sources[0:1]

# For each source in source list
for source,n in zip(sources,range(1,len(sources)+1)):

    # Unpack source
    obsid,srcid,ra,dec,srcnh = source

    try:
        # Print current source
        print('\n\n\nSOURCE %d/%d\n\n\n' % (n,len(sources)))

        # Set model
        model = choose_model(source)
        print('\n\n\nMODEL: %s\n\n\n' % model)

        # Fit sources
        fitting(source,model,bounds=[0.3,10],silent=True)

        # Log stats
        log_info(source,model)

        # Print runtime
        print('\nRUNTIME FOR SOURCE %s-%s: %d sec\n' % (obsid,srcid,time.time()-start_time))


    except:
        print('\n\n\nERROR in SOURCE %s_%s' % (obsid,srcid))

# Print runtime
print('\n\n\nTOTAL RUNTIME: %d\n\n\n' % (time.time()-start_time))
