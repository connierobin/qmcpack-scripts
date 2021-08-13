#! /usr/bin/env python

import optparse
import os
import sys
import copy

import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import math

# Colors to plot with. From ColorBrewer, at colorbrewer2.org
color_list = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3','#ff7f00', '#ffff33', '#a65628', '#f781bf', '#999999']

# Hard coded maximum limit of the x axis
maxIter = 100           # max value of x axis

# Arrays to hold the input info given by the user
qmc_files = []  # array of filenames to read through
plotting_choices = []   # array of arrays of numbers indicating which quantities to plot for each given file
data_labels = []        # a list of the labels to use for the plot

# Strings to search for in the output files
# Different data is printed at each iteration when using the Linear Method or another method
# If certain data is not printed for a method, a nonsense string is used as the indicator for
# the data to ensure that there are no hits for it.
LM_str_dict = {'energy': '  le_mean =  ',
               'uncertainty_of_variance': 'uncertainty =',  # std err of variance
               'variance': 'le_variance = ',
               'std_err': '  stat err =  ',
               'target': '  target function =  ',
               'target_err': '  target stat err = ',
               'std_dev': '  std dev =  ',
               'grad_norm': 'not printed in Linear Method',
               'method_identifier': 'Inside LM engine\'s get_param',
               'data_start': 'After engine_checkConfig',
               'data_end': 'Solving the linear method equations'}

#TODO:Connie look at computeFinalizationUncertainties() in DescentEngine.cpp to see if the uncertainty
# of variance is in fact printed in the descent methods. 
descent_str_dict = {'energy': 'Energy Average:',
                    'uncertainty_of_variance': 'not printed outside Linear Method',
                    'variance': 'Energy Variance:',
                    'std_err': 'Energy Standard Error:',
                    'target': 'Target Function Average:',
                    'target_err': 'Target Function Error:',
                    'std_dev': 'Energy Standard Deviation:',
                    'grad_norm': 'Norm of Gradient Vector:',    # grad norm of target function, not of energy. 
                    'method_identifier': 'Omega from input file',
                    'data_start': 'Before engine_checkConfigurations',
                    'data_end': 'After engine_checkConfig'}

target_dict = {'excited': 'Target: excited state',
               'excited_closest': 'Target: closest excited state',
               'gvp': 'Target: generalized variational principle',
               'ground': 'Target: ground state'}

# TODO:Connie change from curly braces to a regular array
# Arrays to hold the lists of data from each file. 
# If a file does not have this information, the entry will be a blank list for the
# element corresponding to that file. 
energies = {}
uncertainties_of_variances = {}
variances = {}
std_err = {}
target_fn = {}
target_std_err = {}
std_dev = {}
grad_norms = {}

# Values for the average of the last 10 iterations for various values. 
# Will only be correct for the current file being read. 
avg10E = 0.0
avg10std_err = 0.0
avg10target_fn = 0.0
avg10target_std_err = 0.0
avg10std_dev = 0.0
avg10variances = 0.0
avg10uncertainty_of_variances = 0.0
avg10grad_norms = 0.0

# Arrays to hold the lists of data to plot. 
# If data does not exist or is not supposed to be plotted for a given file and data piece,
# the corresponding list will be an empty list. 
energy_to_plot = []
uncertainties_of_variances_to_plot = []
variance_to_plot = []
std_err_to_plot = []
target_fn_to_plot = []
target_std_err_to_plot = []
std_dev_to_plot = []
grad_norms_to_plot = []

# Holds all of the avg10 data for all the files. Each list holds the data in the order of the files. 
# Data for a given data piece for a given file is an average from the last 10 iterations. 
energy_avg10_list = []
uncertainties_of_variances_avg10_list = []
variance_avg10_list = []
std_err_avg10_list = []
target_fn_avg10_list = []
target_std_err_avg10_list = []
std_dev_avg10_list = []
grad_norms_avg10_list = []

# Reads the options provided from the command line. 
def parse_options(args):
    '''Parse arguments from the command line. Expecting 1 input file name.'''

    parser = optparse.OptionParser(usage=__doc__)
    (options, filename) = parser.parse_args(args)

    if len(filename) == 0:
        parser.print_help()
        sys.exit(1)

    return (options, filename)

# Reads the input file to gather the names of the output files to read, 
# the data labels to use, and which sets of data to plot. 
def read_input_file(input_file):
    f = open(input_file, 'r')
    input_lines = f.readlines()

    for line in input_lines:
        words = line.split()
        qmc_files.append(words[0])
        data_labels.append(words[1])
        plotting_choices.append(words[2:])

# Pulls data from the given file and populates all the relevant global variables.
# Namely:
# - LinearMethodUsed or DescentMethodUsed
# - lists: energies, variances, uncertainties_of_variances, std_err, target_fn, target_std_err, std_dev, grad_norms
def extract_data(qmc_file):

    # Need to identify based on the file whether it's using LM or just descent methods
    LinearMethodUsed = False
    DescentMethodUsed = False

    # Try to identify based on the file which type of variational principle is being optimized
    TargetExcitedUsed = False
    TargetExcitedClosestUsed = False
    TargetGVPUsed = False
    TargetGroundUsed = False

    # Set up variables for looping through each line of the file
    have_data = False
    iteration = 0
    f = open(qmc_file)

    for line in f:
        # Check if the identifier line for either is here, to identify which method this file uses
        if not (LinearMethodUsed or DescentMethodUsed):
            if LM_str_dict['method_identifier'] in line:
                LinearMethodUsed = True
                str_dict = LM_str_dict
                print("The following values are reported from the Linear Method")
            elif descent_str_dict['method_identifier'] in line:
                DescentMethodUsed = True
                str_dict = descent_str_dict
                print("The following values are reported from a descent method")
            else:
                continue

        # Check if the target identifier line is used, to identify which variational principle is targeted
        if not (TargetExcitedUsed or TargetExcitedClosestUsed or TargetGVPUsed or TargetGroundUsed):
            if target_dict['excited'] in line:
                TargetExcitedUsed = True
            elif target_dict['excited_closest'] in line:
                TargetExcitedClosestUsed = True
            elif target_dict['gvp'] in line:
                TargetGVPUsed = True
            elif target_dict['ground'] in line:
                TargetGroundUsed = True

        if str_dict['data_start'] in line:
            have_data = True

        if have_data:
            if str_dict['energy'] in line:
                values = line.split()
                energies[iteration] = float(values[2])
            elif str_dict['uncertainty_of_variance'] in line:
                values = line.split()
                uncertainties_of_variances[iteration] = float(values[2])
            elif str_dict['variance'] in line:
                values = line.split()
                variances[iteration] = float(values[2])
            elif str_dict['std_err'] in line:
                values = line.split()
                std_err[iteration] = float(values[3])
            elif str_dict['target'] in line:
                values = line.split()
                if 'N/A' in values[3]:
                    target_fn[iteration] = 0.0
                else:
                    target_fn[iteration] = float(values[3])
            elif str_dict['target_err'] in line:
                values = line.split()
                target_std_err[iteration] = float(values[4])
            elif str_dict['std_dev'] in line:
                values = line.split()
                std_dev[iteration] = float(values[3])
            elif str_dict['grad_norm'] in line:
                grad_norms[iteration] = float(values[3])
            elif str_dict['data_end'] in line:
                # TODO:Connie modify this so that instead of "iteration" it's "number of samples"
                # so that LM results and regular descent results are easier to compare
                have_data = False
                iteration += 1

        if str_dict['data_start'] in line:
            have_data = True

    f.close()

# Print all the data extracted from the data files, as well as the average of the
# last 10 numbers for each quantity to get a finalized value.
def print_data():
    # Print the header.
    sys.stdout.write(' #1._Iteration')
    if len(energies) > 0:
        column_text = '2._Energy'
        sys.stdout.write('%22s' % column_text)
    if len(std_err) > 0:
        column_text = '3._Error'
        sys.stdout.write('%22s' % column_text)
    if len(target_fn) > 0:
        column_text = '4._Target_function'
        sys.stdout.write('%22s' % column_text)
    if len(target_std_err) > 0:
        column_text = '5._Target_error'
        sys.stdout.write('%22s' % column_text)
    if len(std_dev) > 0:
        column_text = '6._Standard_deviation'
        sys.stdout.write('%22s' % column_text)
    if len(variances) > 0:
        column_text = '7._Crude:_Variance'
        sys.stdout.write('%22s' % column_text)
    if len(uncertainties_of_variances) > 0:
        column_text = '8._Crude:_uncertainty_of_Variance'
        sys.stdout.write('%22s' % column_text)
    if len(grad_norms) > 0:
        column_text = '9._Grad_norms'
        sys.stdout.write('%22s' % column_text)
    sys.stdout.write('\n')

    # Extract and print information for each iteration.
    for iter in energies:
        sys.stdout.write('     %9d' % iter)
        if len(energies) > 0:
            sys.stdout.write('   %19.12e' % energies[iter])
        if len(std_err) > 0:
            sys.stdout.write('   %19.12e' % std_err[iter])
        if len(target_fn) > 0:
            sys.stdout.write('   %19.12e' % target_fn[iter])
        if len(target_std_err) > 0:
            sys.stdout.write('   %19.12e' % target_std_err[iter])
        if len(std_dev) > 0:
            sys.stdout.write('   %19.12e' % std_dev[iter])
        if len(variances) > 0:
            sys.stdout.write('   %19.12e' % variances[iter])
        if len(uncertainties_of_variances) > 0:
            sys.stdout.write('   %19.12e' % uncertainties_of_variances[iter])
        if len(grad_norms) > 0:
            sys.stdout.write('   %19.12e' % grad_norms[iter])
        sys.stdout.write('\n')

    # I typically take averages over the last 10 linear method iterations for the values I report
    # as the result of optimization
    # The same idea is applicable to other descent methods

    numIter1 = len(energies)
    # nonlocal maxIter
    # maxIter = numIter1      #TODO:Connie change this to actually be the max samples instead of number of iterations

    print("These values are averaged over the last 10 iterations to come to a final value")

    if len(energies) > 0:
        energies_list = list(energies.values())
        print("Energy")
        avg10E = np.mean(energies_list[numIter1-11:numIter1-1])
        print(avg10E)
        energy_avg10_list.append(avg10E)

    if len(std_err) > 0:
        std_err_list = list(std_err.values())
        print("Energy Uncertainty")
        avg10std_err = np.mean(std_err_list[numIter1-11:numIter1-1])
        print(avg10std_err)
        std_err_avg10_list.append(avg10std_err)

    if len(target_fn) > 0:
        target_fn_list = list(target_fn.values())
        print("Target Function")
        avg10target_fn = np.mean(target_fn_list[numIter1-11:numIter1-1])
        print(avg10target_fn)
        target_fn_avg10_list.append(avg10target_fn)

    if len(target_std_err) > 0:
        target_std_err_list = list(target_std_err.values())
        print("Target Function Uncertainty")
        avg10target_std_err = np.mean(target_std_err_list[numIter1-11:numIter1-1])
        print(avg10target_std_err)
        target_std_err_avg10_list.append(avg10target_std_err)

    if len(std_dev) > 0:
        std_dev_list = list(std_dev.values())
        print("Standard Deviation")
        avg10std_dev = np.mean(std_dev_list[numIter1-11:numIter1-1])
        print(avg10std_dev)
        std_dev_avg10_list.append(avg10std_dev)

    if len(variances) > 0:
        variances_list = list(variances.values())
        print("Variance")
        avg10variances = np.mean(variances_list[numIter1-11:numIter1-1])
        print(avg10variances)
        variance_avg10_list.append(avg10variances)

    if len(uncertainties_of_variances) > 0:
        uncertainties_of_variances_list = list(uncertainties_of_variances.values())
        print("Variance Uncertainty")
        avg10uncertainty_of_variances = np.mean(uncertainties_of_variances_list[numIter1-11:numIter1-1])
        print(avg10uncertainty_of_variances)
        uncertainties_of_variances_avg10_list.append(avg10uncertainty_of_variances)

    if len(grad_norms) > 0:
        grad_norms_list = list(grad_norms.values())
        print("Gradient Norm")
        avg10grad_norms = np.mean(grad_norms_list[numIter1-11:numIter1-1])
        print(avg10grad_norms)
        grad_norms_avg10_list.append(avg10grad_norms)

# Save any data that is chosen to be plotted into data lists so that it can
# be plotted after all the files are processed. 
# Save the list of data points to plot, as well as the "avg10" value, which is the average
# of the last 10 data points. 
def collect_data_for_plots(plotting_choice):
    # Check that enough plotting choices are given, and provide instructions if it is not
    if len(plotting_choice) < 8:
        print("Please specify 8 plotting choices in the input file. ")
        print("After the file name, put a data label for that file, then put a series of 1's and 0's with spaces between corresponding to whether you would like to plot each of the following things: ")
        print("1. Energy")
        print("2. Energy Standard Error")
        print("3. Target Function")
        print("4. Target Function Standard Error")
        print("5. Standard Deviation")
        print("6. Variance")
        print("7. Variance Uncertainty")
        print("8. Gradient Norm")
        print("For example:")
        print("exampleFile.out myData1 1 0 1 0 0 0 0 0")
        print("This would include the data from exampleFile.out on the Energy plot and the Target Function plot, and label it myData1. ")

    # Collect the Energy data
    if int(plotting_choice[0]) == 1:
        energy_to_plot.append(copy.deepcopy(energies))
        subplots_to_show[0] = 1
    else:
        energy_to_plot.append([])
    
    # Collect the Energy Uncertainty data
    if int(plotting_choice[1]) == 1:
        std_err_to_plot.append(copy.deepcopy(std_err))
        subplots_to_show[1] = 1
    else:
        std_err_to_plot.append([])

    # Collect the Target Function data
    if int(plotting_choice[2]) == 1:
        target_fn_to_plot.append(copy.deepcopy(target_fn))
        subplots_to_show[2] = 1
    else:
        target_fn_to_plot.append([])

    # Collet the Target Function Uncertainty data
    if int(plotting_choice[3]) == 1:
        target_std_err_to_plot.append(copy.deepcopy(target_std_err))
        subplots_to_show[3] = 1
    else:
        target_std_err_to_plot.append([])

    # Collect the Standard Deviation data
    if int(plotting_choice[4]) == 1:
        std_dev_to_plot.append(copy.deepcopy(std_dev))
        subplots_to_show[4] = 1
    else:
        std_dev_to_plot.append([])

    # Collect the Variance data
    if int(plotting_choice[5]) == 1:
        variance_to_plot.append(copy.deepcopy(variances))
        subplots_to_show[5] = 1
    else:
        variance_to_plot.append([])

    # Collect the Variance Uncertainty data
    if int(plotting_choice[6]) == 1:
        uncertainties_of_variances_to_plot.append(copy.deepcopy(uncertainties_of_variances))
        subplots_to_show[6] = 1
    else:
        uncertainties_of_variances_to_plot.append([])

    # Collect the Gradient Norm data
    if int(plotting_choice[7]) == 1:
        grad_norms_to_plot.append(copy.deepcopy(grad_norms))
        subplots_to_show[7] = 1
    else:
        grad_norms_to_plot.append([])

# Put together the plot using matplotlib, and generate all the needed subplots
# for each quantity that has been scraped. 
def plot_data():
     # TODO:Connie Use the dict names for maximum generality?

    # set up plot
    fig = plt.figure()

    # Add ENERGY subplot
    add_subplot(fig, energy_to_plot, energy_avg10_list, 'Energy', 'Energy (Ha)', 'Energy')

    # Add STANDARD ERROR subplot
    add_subplot(fig, std_err_to_plot, std_err_avg10_list, 'Standard Error', 'Std Err (Ha)', 'Std Err')

    # Add TARGET FUNCTION subplot
    add_subplot(fig, target_fn_to_plot, target_fn_avg10_list, 'Target Function', 'Target Fn Value', 'Targ Fn')

    # Add TARGET FUNCTION STANDARD ERROR subplot
    add_subplot(fig, target_std_err_to_plot, target_std_err_avg10_list, 'Target Function Standard Error', 'Target Std Err', 'Targ Std Err')

    #TODO:Connie be more clear in the plot title what this is the std dev of. I think it's the energy?
    # Add STANDARD DEVIATION subplot
    add_subplot(fig, std_dev_to_plot, std_dev_avg10_list, 'Standard Deviation', 'Std Dev', 'Std Dev')

    # Add VARIANCE subplot
    add_subplot(fig, variance_to_plot, variance_avg10_list, 'Variance', 'Variance', 'Var')

    #TODO:Connie change this to stay std dev or std err instead of "uncertainty"
    # Add TARGET FUNCTION STANDARD ERROR subplot
    add_subplot(fig, uncertainties_of_variances_to_plot, uncertainties_of_variances_avg10_list, 'Variance Uncertainty', 'Var Uncertainty', 'Var Unc')

    #TODO:Connie make it clear what this is the gradient norm of! I think it's the energy?
    # Add GRADIENT NORM subplot
    add_subplot(fig, grad_norms_to_plot, grad_norms_avg10_list, 'Gradient Norm', 'Grad Norm', 'Grad Norm')


    # Show only "outer" axis labels
    for ax in fig.get_axes():
        ax.label_outer()

    # Plot all subplots!
    plt.show()

# TODO:Connie change the x axis to number of samples instead?
# Adds a subplot to the matplotlib figure provided. 
# Plots the values listed in values_to_plot (which is a list of lists of data, some of which may be empty)
# Draws a box to display the values from values_avg10_list
# Places the plot title given from plot_title
# Places the y axis label given from yaxis_label
# Labels the data in a legend using the results_label provided
def add_subplot(fig, values_to_plot, values_avg10_list, plot_title, yaxis_label, results_label):
    added_subplot = False
    results_lines = [] # temporary string to hold the finalized values to be displayed on each plot

    # Add the Std Err data to the subplot
    for i in range(len(values_to_plot)):
        data_list = values_to_plot[i]
        if len(data_list) > 0:
            if added_subplot == False:
                # make sure this only gets run once for subplot setup
                added_subplot = True
                # alter subplot geometry
                n = len(fig.axes)
                for j in range(n):
                    fig.axes[j].change_geometry(n+1, 1, j+1)
                # create subplot
                ax = fig.add_subplot(n+1, 1, n+1)
            # Add data to subplot
            ax.plot(range(1, len(data_list)+1), list(data_list.values()), marker='o', color=color_list[i % len(color_list)],label='%s: %s' % (results_label, data_labels[i]))
            results_lines.append('Avg10 %s: $%.5f$  %s' % (results_label, values_avg10_list[i], data_labels[i]))

    if added_subplot == True:
        # White background for the plot
        props = dict(boxstyle='square', facecolor='white')

        # Axes and title settings
        ax.set_title(plot_title)
        ax.set_xlabel('Iteration Number', horizontalalignment='right')
        ax.set_ylabel(yaxis_label)
        # ax.set_ylim(linE-.01,linE+.01)   # Often helpful to zoom in on converged part of optimization
        ax.get_yaxis().get_major_formatter().set_useOffset(False)

        # Add vertical grid lines
        ax = plt.gca()
        ax.grid(which='major', axis='x', linestyle='--')
        ax.xaxis.set_minor_locator(MultipleLocator(10))
        ax.grid(which='minor', axis='x', linestyle=':')

        # Place the legend
        ax.legend(loc='upper right', bbox_to_anchor=(.95, 1.0))

        # Make a text box for the results string, showing the average of the last 10 energies
        results_str = '\n'.join(results_lines)
        ax.text(.1, .95, results_str, transform=ax.transAxes, horizontalalignment='left', verticalalignment='top', bbox=props)
        
        # Adjust the position and size of the subplot
        box = ax.get_position()
        ax.set_position([box.x0, box.y0+box.height*.15, box.width, box.height*.85])


if __name__ == '__main__':
    (options, input_files) = parse_options(sys.argv[1:])
    read_input_file(input_files[0])
    #energies, uncertainties_of_variances, variances, std_err, target_fn, target_std_err, std_dev, grad_norms = extract_data(data_files)
    for i in range(len(qmc_files)):
        qmc_file = str(qmc_files[i])
        plotting_choice = plotting_choices[i]
        extract_data(qmc_file)
        print_data()
        collect_data_for_plots(plotting_choice)
    plot_data()
