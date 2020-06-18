#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# FILE: acab_funcs.py                                                       #
#                                                                           #
# AUTHOR/UPDATED: Mason V. Tea / 6.15.20                                    #
#                                                                           #
# PURPOSE: Functions for Automated Characterization of Accreting Binaries.  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Import libraries
import os
from sherpa.astro.ui import *
import matplotlib.pylab as plt
import numpy as np
import logging
from itertools import combinations

# Matplotlib style sheet
plt.style.use('bmh')






#--------------------------------------------------------------------------------------------------------------------#
# FUNCTION:     log_info()                                                                                           #
#                                                                                                                    #
# DESCRIPTION:  Logs statistics for a given source with the fit most recently run.                                   #
#                                                                                                                    #
# VARIABLES:    source -- List of five strings; [<obsid>, <source number>, <RA>, <DEC>, <Goddard nH>]                #
#               model -- String; '(<absorption models>)*(emission models)'; suffixes .abs<n> & .emis<n>, n = 1,2,... #
#--------------------------------------------------------------------------------------------------------------------#

def log_info(source,model):

    # Unpack source array
    obsid,srcid,ra,dec,srcnh = source

    # Set logfile for given source
    logfile = 'plots/logs/%s_%s-fit.log' % (obsid,srcid)

    # This process sometimes raises an error...
    try:
        # Calculate characteristic values
        net_src, net_bkg = calc_data_sum(id=1), calc_data_sum(bkg_id=1)
        rate_src, rate_bkg = calc_data_sum()/get_exposure(id=1), calc_data_sum(bkg_id=1)/get_exposure(bkg_id=1)
        pflux, eflux = calc_photon_flux(), calc_energy_flux()

        # Write characteristic values, statistics to logfile for given source
        with open(logfile,'w') as f:

            f.write('Source: %s-%s' % (obsid,srcid))
            f.write('\n(RA,Dec) = (%s,%s)' % (ra,dec))
            f.write('\nModel: %s' % model)
            f.write('\nMilky Way nH = %s\n' % srcnh)

            f.write('\nNet counts (src): ' % net_src)
            f.write(str(net_src))

            f.write('\nCount rate (src): ' % rate_src)
            f.write(str(rate_src))

            f.write('\nNet counts (bkg):' % net_bkg)
            f.write(str(net_bkg))

            f.write('\nCount rate (bkg): ' % rate_bkg)
            f.write(str(rate_bkg))

            f.write('\nPhoton flux: ' % pflux)
            f.write(str(pflux))

            f.write('\nEnergy flux: ')
            f.write(str(eflux))
            f.write('\n\n')

            covariance()
            f.write(str(get_covar_results()))

            f.write(str(get_stat_info()[0]))

        f.close()

    # Log *only* error message if error is raised
    except:
        f = open(logfile,'w')
        f.write('Error in calculations. Check PHA files.')
        f.close()






#--------------------------------------------------------------------------------------------------------------------#
# FUNCTION:     fitting_silent()                                                                                     #
#                                                                                                                    #
# DESCRIPTION:  Identical to fitting(), with print statements removed.                                               #
#                                                                                                                    #
# VARIABLES:    source -- List of five strings; [<obsid>, <source number>, <RA>, <DEC>, <Goddard nH>]                #
#               model -- String; '(<absorption models>)*(emission models)'; suffixes .abs<n> & .emis<n>, n = 1,2,... #
#               binning -- Integer; counts per bin                                                                   #
#               bounds -- List of two floats; [<lower bound>, <upper bound>]                                         #
#               silent -- Boolean; allow or disallow fitting() to print statements to terminal                       #
#--------------------------------------------------------------------------------------------------------------------#

def fitting_silent(source, model, binning='', bounds=''):

        # Unpack source array to define OBSID, source ID, column density
        obsid,srcid,srcnh = source[0],source[1],source[4]

        #------------------------#
        #     INITIALIZATION     #
        #------------------------#

        # Load grouped spectrum
        load_pha('%s/extracted_spectra_%s_%s_grp.pi' % (obsid,obsid,srcid))

        #------------------------#
        #      RAW SPECTRUM      #
        #------------------------#

        # Save plot of raw spectrum in log space
        set_xlog()
        set_ylog()
        plot_data()
        plt.savefig('plots/raw/%s_%s-rawspec.pdf' % (obsid,srcid))
        plt.close()

        #-----------------------#
        #     ENERGY BOUNDS     #
        #-----------------------#

        # If bounds are not provided with function call, prompt for user input
        if bounds == '':

            # Open raw spectrum
            os.system('open plots/raw/%s_%s-rawspec.pdf' % (obsid,srcid))

            # Error if user hits 'enter,' so pass if error is raised (no energy bounds attributed in this case)
            try:
                lo,hi = input('\nfitting(): Enter energy bounds, separated by space (Enter to skip): ').split()
                notice(float(lo),float(hi))
                plot_data()
                plt.savefig('plots/raw/%s_%s-rawspec_bounded.pdf' % (obsid,srcid))
                plt.close()
            except:
                pass

        # If bounds are given in function call, use those
        else:
            notice(bounds[0],bounds[1])
            plot_data()
            plt.savefig('plots/raw/%s_%s-rawspec_bounded.pdf' % (obsid,srcid))
            plt.close()

        #-----------------------#
        #        BINNING        #
        #-----------------------#

        # If no binning is provided in function call, use default binning
        if binning == '':
            pass

        # If binning is provided in function call, use the provided counts/bin value
        else:
            group_counts(binning)

        #----------------------#
        #         MODEL        #
        #----------------------#

        # Subtract background
        subtract()

        # Set model
        set_source(model)

        # Set statistic
        set_stat("chi2datavar")

        # Set hydrogen column density, abundance & cross-section for xstbabs
        abs1.nH = float(srcnh)/1.0E22
        freeze(abs1.nH)
        set_xsabund('wilm')
        set_xsxsect('vern')

        #----------------------#
        #        FITTING       #
        #----------------------#

        # Apply fit
        fit()

        # Plot, save and open fit
        plot_fit_delchi()
        plt.savefig('plots/fits/%s_%s_fit.pdf' % (obsid,srcid))
        plt.close()







#--------------------------------------------------------------------------------------------------------------------#
# FUNCTION:     fitting()                                                                                            #
#                                                                                                                    #
# DESCRIPTION:  Generates fit of grouped Chandra spectrum using SHERPA.                                              #
#                                                                                                                    #
# VARIABLES:    source -- List of five strings; [<obsid>, <source number>, <RA>, <DEC>, <Goddard nH>]                #
#               model -- String; '(<absorption models>)*(emission models)'; suffixes .abs<n> & .emis<n>, n = 1,2,... #
#               binning -- Integer; counts per bin                                                                   #
#               bounds -- List of two floats; [<lower bound>, <upper bound>]                                         #
#               silent -- Boolean; allow or disallow fitting() to print statements to terminal                       #
#--------------------------------------------------------------------------------------------------------------------#

def fitting(source, model, binning='', bounds='', silent=False):

        # Unpack source array to define OBSID, source ID, column density
        obsid,srcid,srcnh = source[0],source[1],source[4]

        # Run with print statements
        if silent == False:

            #------------------------#
            #     INITIALIZATION     #
            #------------------------#

            # Load grouped spectrum
            load_pha('%s/extracted_spectra_%s_%s_grp.pi' % (obsid,obsid,srcid))

            #------------------------#
            #      RAW SPECTRUM      #
            #------------------------#

            # Save plot of raw spectrum in log space
            print('\nfitting(): Plotting raw spectrum...')
            set_xlog()
            set_ylog()
            plot_data()
            plt.savefig('plots/raw/%s_%s-rawspec.pdf' % (obsid,srcid))
            plt.close()
            print('\nfitting(): Saved as %s_%s-rawspec.pdf.' % (obsid,srcid))

            #-----------------------#
            #     ENERGY BOUNDS     #
            #-----------------------#

            # If bounds are not provided with function call, choose default inputß
            if bounds == '':
                bounds = [0.3,10]
            else:
                pass

            notice(bounds[0],bounds[1])
            plot_data()
            print('\nfitting(): Energy bounds set to [' + str(bounds[0]) + ' ,' + str(bounds[1]) + '].')
            print('\nfitting(): Saving bounded spectrum as %s_%s-rawspec_bounded.pdf' % (obsid,srcid))
            plt.savefig('plots/raw/%s_%s-rawspec_bounded.pdf' % (obsid,srcid))
            plt.close()

            #-----------------------#
            #        BINNING        #
            #-----------------------#

            # If no binning is provided in function call, use default binning
            if binning == '':
                pass

            # If binning is provided in function call, use the provided counts/bin value
            else:
                group_counts(binning)
                print('\nfitting(): Binning set to %d counts per bin.' % binning)

            #----------------------#
            #         MODEL        #
            #----------------------#

            # Subtract background
            subtract()

            # Set model
            print('\nfitting(): Fitting spectrum with model %s...' % model)
            set_source(model)

            # Set statistic
            set_stat("chi2datavar")

            # Set hydrogen column density, abundance & cross-section for xstbabs
            abs1.nH = float(srcnh)/1.0E22
            freeze(abs1.nH)
            set_xsabund('wilm')
            set_xsxsect('vern')

            #----------------------#
            #        FITTING       #
            #----------------------#

            # Apply fit
            print('\nfitting(): Fitting...')
            fit()

            # Plot, save and open fit
            print('\nfitting(): Plotting fit...\n')
            plot_fit_delchi()
            print('fitting(): Printing stats...\n')
            print('fitting(): Saving figure as %s_%s_fit.pdf' % (obsid,srcid))
            plt.savefig('plots/fits/%s_%s_fit.pdf' % (obsid,srcid))
            plt.close()

        # Run without print statements
        else:
            fitting_silent(source,model,binning,bounds)






#--------------------------------------------------------------------------------------------------------------------#
# FUNCTION:     fitting_fast()                                                                                       #
#                                                                                                                    #
# DESCRIPTION:  Identical to fitting_silent(), with plot statements removed.                                         #
#                                                                                                                    #
# VARIABLES:    source -- List of five strings; [<obsid>, <source number>, <RA>, <DEC>, <Goddard nH>]                #
#               model -- String; '(<absorption models>)*(emission models)'; suffixes .abs<n> & .emis<n>, n = 1,2,... #
#               binning -- Integer; counts per bin                                                                   #
#               bounds -- List of two floats; [<lower bound>, <upper bound>]                                         #
#--------------------------------------------------------------------------------------------------------------------#

def fitting_fast(source, model, binning='', bounds=''):

        # Unpack source array to define OBSID, source ID, column density
        obsid,srcid,srcnh = source[0],source[1],source[4]

        #------------------------#
        #     INITIALIZATION     #
        #------------------------#

        # Load grouped spectrum
        load_pha('%s/extracted_spectra_%s_%s_grp.pi' % (obsid,obsid,srcid))

        #-----------------------#
        #     ENERGY BOUNDS     #
        #-----------------------#

        # If bounds are not provided with function call, use default bounds
        if bounds == '':
            bounds = [0.3,10]
        else:
            pass

        # Set bounds
        notice(bounds[0],bounds[1])

        #-----------------------#
        #        BINNING        #
        #-----------------------#

        # If no binning is provided in function call, use default binning
        if binning == '':
            pass

        # If binning is provided in function call, use the provided counts/bin value
        else:
            group_counts(binning)

        #----------------------#
        #         MODEL        #
        #----------------------#

        # Subtract background
        subtract()

        # Set model
        set_source(model)

        # Set statistic
        set_stat("chi2datavar")

        # Set hydrogen column density, abundance & cross-section for xstbabs
        abs1.nH = float(srcnh)/1.0E22
        freeze(abs1.nH)
        set_xsabund('wilm')
        set_xsxsect('vern')

        #----------------------#
        #        FITTING       #
        #----------------------#

        # Apply fit
        fit()

        # Plot, save and open fit
        plot_fit_delchi()
        plt.savefig('plots/fits/%s_%s_fit.pdf' % (obsid,srcid))
        plt.close()







#--------------------------------------------------------------------------------------------------------------------#
# FUNCTION:     ftest()                                                                                              #
#                                                                                                                    #
# DESCRIPTION:  Perform statistical comparison (f-test) between two models, return best model.                       #
#                                                                                                                    #
# VARIABLES:    source -- List of five strings; [<obsid>, <source number>, <RA>, <DEC>, <Goddard nH>]                #
#               model1 -- String; '<model>'                                                                          #
#               model2 -- String; '<model>'                                                                          #
#               binning -- Integer; counts per bin                                                                   #
#               bounds -- List of two floats; [<lower bound>, <upper bound>]                                         #
#               good -- Boolean; determines whether good or bad model is returned                                    #
#--------------------------------------------------------------------------------------------------------------------#



def ftest(source,model1,model2,binning='',bounds='',good=True):

    # Fit simple model, get reduced chi-squared and degrees of freedom
    fitting_fast(source, model1, binning, bounds)
    chi1, dof1 = get_stat_info()[0].rstat, get_stat_info()[0].dof

    # Fit simple model, get reduced chi-squared and degrees of freedom
    fitting_fast(source, model2, binning, bounds)
    chi2, dof2 = get_stat_info()[0].rstat, get_stat_info()[0].dof

    # Decide which model is more complex
    if dof1 > dof2:
        simpMod, compMod = model1, model2
        simpChi, compChi = chi1, chi2
        simpDOF, compDOF = dof1, dof2
    elif dof2 > dof1:
        simpMod, compMod = model2, model1
        simpChi, compChi = chi2, chi1
        simpDOF, compDOF = dof2, dof1
    # If DoF vals are the same, calc_ftest doesn't work; arbitrarily choose the model with the lower red. chi^2 in this case
    elif chi1 < chi2:
        with open('tournament.log','a') as f: f.write('\n\n%s BEATS %s\n\n' % (model1,model2))
        f.close()
        return(model1)
    # If red. chi^2 are identical, arbitrarily choose model2
    elif chi2 <= chi1:
        with open('tournament.log','a') as f: f.write('\n\n%s BEATS %s (TOSS UP)\n\n' % (model2,model1))
        f.close()
        return(model2)

    # Calculate p-value, select model (95% certainty)
    p = calc_ftest(simpDOF, simpChi, compDOF, compChi)

    # Return model, log results
    if p > 0.05:
        if good==True:
            with open('tournament.log','a') as f: f.write('\n\n%s BEATS %s\n\n' % (simpMod,compMod))
            f.close()
            return(simpMod)
        else:
            with open('tournament.log','a') as f: f.write('\n\n%s BEATS %s\n\n' % (compMod,simpMod))
            f.close()
            return(compMod)
    else:
        if good==True:
            with open('tournament.log','a') as f: f.write('\n\n%s BEATS %s\n\n' % (compMod,simpMod))
            f.close()
            return(compMod)
        else:
            with open('tournament.log','a') as f: f.write('\n\n%s BEATS %s\n\n' % (simpMod,compMod))
            f.close()
            return(simpMod)







#--------------------------------------------------------------------------------------------------------------------#
# FUNCTION:     choose_model()                                                                                       #
#                                                                                                                    #
# DESCRIPTION:  Perform statistical comparison (f-test) between all available models, return best model.             #
#                                                                                                                    #
# VARIABLES:    source -- List of five strings; [<obsid>, <source number>, <RA>, <DEC>, <Goddard nH>]                #
#               binning -- Integer; counts per bin                                                                   #
#               bounds -- List of two floats; [<lower bound>, <upper bound>]                                         #
#--------------------------------------------------------------------------------------------------------------------#

def choose_model(source,binning='',bounds=''):

    #--------------------------------#
    #        MODEL DEFINITIONS       #
    #--------------------------------#

    # Define default base model (always two-component absorption from MW & source)
    base = '(xstbabs.abs1+xstbabs.abs2)'

    # Define allowed emission models
    diskModels = ['xsdiskbb.didkbb','xsdiskpn.diskpn','xsapec.apec']
    coronaModels = ['powlaw1d.powlaw','xscompbb.compbb']
    gasModels = ['xsmekal.mekal']

    # Create an array of all models; gas, corona cannot exist without disk & base model alone not sufficient
    allModels = []
    for disk in diskModels:
        for corona in coronaModels:
            for gas in gasModels:
                # allModels.append(base + ('*(%s)' % (gas)))
                # allModels.append(base + ('*(%s+%s)' % (gas,corona)))
                allModels.append(base + ('*(%s+%s)' % (gas,disk)))
                allModels.append(base + ('*(%s+%s+%s)' % (gas,corona,disk)))
            # allModels.append(base + ('*(%s)') % corona)
            allModels.append(base + ('*(%s+%s)' % (corona,disk)))
        # allModels.append(base + '*(%s)' % disk)

    #------------------------#
    #        F-TESTING       #
    #------------------------#

    # Loop until 'winner' is found
    winner = 0
    contenders = allModels
    while winner == 0 or len(contenders) > 1:

        # Create list of tuples containing possible models for comparison
        combos = list(combinations(contenders,2))

        # For each combination, run an f-test and remove the loser from the list of possible models
        for combo in combos:
            print('\n\n\n\nLEN(COMBOS) = %d, LEN(CONTENDERS) = %d\n\n\n\n' % (len(combos),len(contenders)))
            badmodel = ftest(source,combo[0],combo[1],binning=binning,bounds=bounds,good=False)
            for n,model in zip(range(len(contenders)),contenders):
                if badmodel == model:
                    contenders.pop(n)
            for m,combo2 in zip(range(len(combos)),combos):
                if badmodel == combo2[0] or badmodel == combo2[1]:
                    combos.pop(m)

        # When one model remains, exit the loop
        if len(contenders) > 1:
            winner = 0
        else:
            # Return best model
            return(contenders[0])

    # Return best model
    return(contenders[0])
