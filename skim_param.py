import os, sys
import matplotlib.pyplot as plt
import numpy as np
import statistics


# Get the norm of a parameter set and find the relative importance
# Of different pieces.  
# Right now hard coded to look for how important the first three
# Determinants in the CI expansion
def norm_param(list_coeff):
    """
    Get the norm of a parameter set and find the relative importance
    Of different pieces.  
    Right now hard coded to look for how important the first three
    Determinants in the CI expansion
    """
    total = sum([coeff*coeff for coeff in list_coeff])
    print("First 3:")
    for i in range(3):
       print(list_coeff[i]*list_coeff[i]/total)
    
# Input is the raw QMC Output file
filename = sys.argv[1]
print(filename)


    
updatesection = False
collecting = False
havetags=False

this_iter_vars = []
all_iters_vars=[]
tags=[]

f=open(filename,'r')
for line in f:
    # Collecting set to true upon finding an <optVariables> tag
    # Stop collecting when we hit the end of the variable section
    if collecting and "</optVariables" not in line:
        this_iter_vars.append(float(line.split()[1]))
        # Also grab the tags for the different variable sets (uu, ud, eX, ci, etc.)
        # Only do this on the first iteration
        if not havetags:
            temp = str(line.split()[0])
            tags.append(temp.split("_")[0])
    elif "<optVariables" in line and not updatesection:
        collecting = True
    elif "</optVariables" in line:
        if len(tags)>0:
            havetags=True
        collecting = False
        if (len(this_iter_vars)>0):
            all_iters_vars.append(this_iter_vars)
        this_iter_vars = []
        if updatesection:
            updatesection = False
    # If using the hybrid method, linear method, or block linear method
    # An intermediate set of parameters is printed, don't collect these
    elif "initial energy" in line or "Updating the guiding" in line:
        updatesection=True
    elif "Applying the update" in line:
        updatesection=False
f.close()

# Transpose the list of list of optimizable parameters, so plotting is simpler
np_vars = (np.array(all_iters_vars)).T

# Each tag repeats many times. Use a list of unique tags to identify the ranges of tags that
# are all the same tag, and store things tag -> range mappings in var_sets. 
unique_tags=list(set(tags))
var_sets={}
one_set=[]
for tag in unique_tags:
    #Get left and right bounds of tag indicies
    left = min([i for i in range(len(tags)) if tags[i]==tag])
    right = max([i for i in range(len(tags)) if tags[i]==tag])+1
    var_sets[tag]=np_vars[left:right]

# Make a separate plot for each tag
for key in var_sets:
    plt.figure(key)
    item=var_sets[key]
    if "CI" in key:
        quick_fix = np.array(item).T
        norm_param(quick_fix[len(quick_fix)-1])
    [plt.plot(np.arange(len(this_var)),this_var) for this_var in item]
    plt.title(key)
plt.show()

