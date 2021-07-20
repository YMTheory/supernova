import sys
sys.path.append('/dybfs2/nEXO/fuys/2020-sensitivity/sensitivity/modules')

import pandas as pd
import histlite as hl
import numpy as np
from matplotlib import pyplot as plt

plt.rcParams.update({'font.size': 14})
plt.rcParams['figure.figsize'] = (7,6)

import nEXOFitWorkspace
import nEXOFitModel
import nEXOFitLikelihood
#1 Creating a workspace
workspace = nEXOFitWorkspace.nEXOFitWorkspace(config='/dybfs2/nEXO/fuys/sensitivity/work/config/TUTORIAL_config.yaml')
#workspace = nEXOFitWorkspace.nEXOFitWorkspace(config='../config/Sensitivity2020_config_customBinningTest.yaml')
workspace.LoadComponentsTableFromFile( 'ComponentsTable_D-005_V2019.h5' )
workspace.CreateGroupedPDFs()

likelihood = nEXOFitLikelihood.nEXOFitLikelihood()
likelihood.AddPDFDataframeToModel( workspace.df_group_pdfs, workspace.histogram_axis_names )
initial_guess = likelihood.GetVariableValues()


# Next, set limits so that none of the PDFs go negative in the fit.
for var in likelihood.model.variable_list:
    if 'Bb0n' in var['Name']:
        likelihood.SetVariableLimits( var['Name'], \
                                  lower_limit = 0., \
                                  upper_limit = 80.)
    else: 
        likelihood.SetVariableLimits( var['Name'], \
                                  lower_limit = 0., \
                                  upper_limit = var['Value']*10.)
        
        
# Fluctuate and set Rn222 constraint
rn222_idx = likelihood.GetVariableIndex('Rn222')
rn222_constraint_val = (np.random.randn()*0.1 + 1)*initial_guess[rn222_idx]
likelihood.SetGaussianConstraintAbsolute(likelihood.model.variable_list[rn222_idx]['Name'],\
                                                         rn222_constraint_val, \
                                                         0.1 * initial_guess[rn222_idx])


# Fix the Co60 parameter, since this PDF is really underconstrained.
likelihood.SetVariableFixStatus('Num_FullTPC_Co60',True)

# Increase the step size for the Bb0n variable
likelihood.SetFractionalMinuitInputError('Num_FullLXeBb0n', 0.01/0.0001)

print('\n\n')
likelihood.PrintVariableList() 

print('\nConstraints:')
for constraint in likelihood.model.constraints:
    print('\t{}'.format(constraint))
print('\n')

num_datasets = 5   # How many datasets to generate for our ensemble
num_hypotheses = 10    # How many different values of N_signal to try
xvals = np.linspace(0,28,num_hypotheses)    # Array of values to be used as N_signal hypotheses
lambdas = np.zeros((num_datasets,num_hypotheses))   # Empty array which will contain the computed test statistics
best_fit_converged = np.ones(num_datasets,dtype=bool)    # Array for flagging failed fits
converged_mask = np.ones((num_datasets,num_hypotheses),dtype=bool)

global_initial_values = np.copy(initial_guess)

import time
start = time.time() # Just for benchmarking

# Set the verbosity. Heads up - this prints out a LOT of stuff
VERBOSE = False

for j in range(0,num_datasets):
    
    print('Running dataset {} at {:3.3} min\n'.format(j,(time.time()-start)/60.))
    
    initial_values = np.copy(global_initial_values)
    
    # Create a new dataset, using the inital input parameters (no signal.)
    likelihood.model.UpdateVariables( initial_values )
    likelihood.model.GenerateModelDistribution()
    likelihood.AddDataset( likelihood.model.GenerateDataset() )
    
    # Fix the Co60 paramter, since it's underconstrained and can lead to failures
    # the fit convergence.
    likelihood.SetAllVariablesFloating()
    likelihood.SetVariableFixStatus('Num_FullTPC_Co60',True)
    
    # Fluctuate and set the Rn222 constraint, (the fluctuation avoids a bias in the
    # pull distributions)
    rn222_idx = likelihood.GetVariableIndex('Rn222')
    rn222_constraint_val = (np.random.randn()*0.1 + 1)*initial_guess[rn222_idx]
    likelihood.SetGaussianConstraintAbsolute(likelihood.model.variable_list[rn222_idx]['Name'],\
                                                         rn222_constraint_val, \
                                                         0.1 * initial_guess[rn222_idx])


    if VERBOSE:
                likelihood.PrintVariableList()
                print('\nConstraints:')
                for constraint in likelihood.constraints:
                    print('\t{}'.format(constraint))
                print('\n')

        
    for i in (range(0,num_hypotheses)):
            
            print('Running hypothesis {} ({:3.3} signal counts)'.format(i,xvals[i]))
            
            signal_idx = likelihood.GetVariableIndex( 'Num_FullLXeBb0n' ) # Get the index of the 0nu variable
            
            initial_values = np.copy(global_initial_values)
            initial_values[signal_idx] = xvals[i]
    
            # Here, we fix the number of signal events to a specific hypothesis.
            likelihood.SetVariableFixStatus('Bb0n',True)
            
            lambda_fit_result = likelihood.ComputeLambdaForPositiveSignal(\
                                                        initial_values=initial_values,\
                                                        signal_name='Bb0n',\
                                                        signal_expectation=xvals[i],\
                                                        fixed_fit_signal_value=xvals[i],\
                                                        print_level=0 )
            
            lambdas[j,i]   = lambda_fit_result['lambda']
            num_iterations = lambda_fit_result['fixed_fit_iterations']
            converged_mask[j][i] = lambda_fit_result['fixed_fit_covar'] # if the fixed fit had a good covariance matrix
            
            if VERBOSE:
                print('\tLambda = {:3.3}'.format(lambdas[j,i]))
                print('\tNumber of iterations: {}'.format(int(num_iterations)))
                print('\tIs valid? {}'.format(likelihood.fitter.get_fmin()['is_valid']))
                print('\tAccurate covariance? {}'.format(likelihood.fitter.get_fmin()['has_accurate_covar']))
                print('\tForced pos def covariance? {}'.format(likelihood.fitter.get_fmin()['has_made_posdef_covar']))
                #likelihood.PrintVariableList()
                print('\n')
        
    print('\t{} fits converged ({:4.4}%)\n\n'.format( np.sum(converged_mask[j]), \
                                    100.*np.sum(converged_mask[j])/len(converged_mask[j])) )
                
print('\nElapsed: {:4.4}s'.format(time.time()-start))
plt.rcParams.update({'font.size': 14})
plt.rcParams['figure.figsize'] = (7,6)

for i in range(0,num_datasets):
    
    # Note: I only plot points for which the fit converged. 
    if best_fit_converged[i]:
        plt.plot(xvals[converged_mask[i]],lambdas[i][converged_mask[i]],'-',label='Run{}'.format(i),linewidth=0.5)
        #plt.plot(xvals[converged_mask[i]],lambdas[i][converged_mask[i]],'-o',label='Run{}'.format(i),markersize=3)
plt.plot(xvals,np.ones(len(xvals))*2.706,'--k')
plt.text(7.,2.9,'90% confidence threshold (Wilks\')',fontsize=14)
plt.axis([0.,28.,0.,10.])
#plt.legend(loc='upper left',fontsize=12,ncol=2,labelspacing=0.1)
plt.ylabel('Test statistic')
plt.xlabel('Number of 0nuBB counts')
plt.savefig('/dybfs2/nEXO/fuys/2020-sensitivity/sensitivity/work/SensitivityPaper2020_scripts/nubb_counts.png')
print('Number in which the best-fit converged: {}'.format(np.sum(best_fit_converged)))

curve_idx = 1 # Choose one of the curves, it doesn't matter which one

plt.figure()
plt.plot(xvals,lambdas[curve_idx],'-o',label='Example toy experiment')
plt.plot(xvals,np.ones(len(xvals))*2.706,'--k')

# Define a mask around the crossing point, and use only fits which converged:
indexes = range(len(lambdas[curve_idx]))
mask = (lambdas[curve_idx]>1.)&\
        (lambdas[curve_idx]<7.)&\
        (converged_mask[curve_idx])&\
        (indexes>np.argmin(lambdas[curve_idx]))

if len(xvals[mask]) > 0:
    # Fit a quadratic function, then plot it on top of the test statistic curve
    p = np.polyfit(xvals[mask],lambdas[curve_idx][mask],2.)
    plt.plot(xvals, p[0]*xvals**2 + p[1]*xvals + p[2],'--',linewidth=2,label='Quadratic fit')

plt.text(9.,2.9,'90% confidence threshold',fontsize=14)
plt.text(10.5,2.3,'(Wilks\' theorem)',fontsize=14)
plt.axis([0.,28.,0.,6.])
plt.legend(loc='upper left',fontsize=11)
plt.ylabel('Likelihood ratio test statistic')
plt.xlabel('Assumed number of 0nuBB counts')
plt.savefig('./0nbb_quadratic_fit')


crossings = np.zeros(num_datasets) - 1.

for i in range(0,num_datasets):
    
        # only fit around the threhsold...
        mask = (lambdas[i]>1.)&(lambdas[i]<7.)&converged_mask[i]
        
        if len(xvals[mask]) > 0:
            p = np.polyfit(xvals[mask],lambdas[i][mask],2.)
            crossings[i] = (-p[1] + np.sqrt( p[1]**2 - 4*(p[0])*(p[2]-2.706) ))/(2*p[0])
        
print(np.mean(crossings),'---mean----')
print(np.median(crossings),'---median----')

fig,ax = plt.subplots(1,1)
hteststats = hl.hist(crossings[crossings>0.],bins=np.linspace(0,30,31))
hl.plot1d(ax,hteststats,color='r')
ax.set_ylabel('Counts ({} trials total)'.format(num_datasets))
ax.set_xlabel('0nuBB counts at 90% confidence limit')
plt.savefig('./example_sensitivity_distribution.png',dpi=400,bbox_inches='tight')

plt.figure()
hcdf = hl.Hist(hteststats.bins,np.cumsum(hteststats.values)/np.sum(hteststats.values))
hl.plot1d(hcdf,label='Normalized CDF')
plt.plot(np.linspace(0.,30.,10),np.ones(10)*0.5,'--k',label='Median')

plt.axis([0.,30.,0.,1.1])
plt.grid()
plt.legend(loc='lower right')

plt.ylabel('Fractional cumulative counts')
plt.xlabel('0nuBB counts at 90% confidence limit')
plt.savefig('./Fractional_cumulative.png',dpi=400,bbox_inches='tight')
