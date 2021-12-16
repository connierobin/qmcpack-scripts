import os, sys
import matplotlib.pyplot as plt
import numpy as np
import statistics


# Get the norm of a list and find the relative importance
# of different pieces.  
# Right now hard coded to look for how important the first three
# values in the list
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
    return [coeff*coeff / total for coeff in list_coeff]
    
# Input is the raw QMC Output file
filename = sys.argv[1]
print(filename)

collecting_eDerivs = False
collecting_lDerivs = False

this_iter_ederivs = []
all_iters_ederivs=[]
this_iter_lderivs = []
all_iters_lderivs=[]

f=open(filename,'r')
for line in f:
    # Collecting set to true upon finding an <eDerivs> tag
    # Stop collecting when we hit the end of the variable section
    if not collecting_eDerivs and not collecting_lDerivs:
        if "<e_derivs" in line:
            collecting_eDerivs = True
        if "<lderivs" in line:
            collecting_lDerivs = True
    else:
        if collecting_eDerivs and "</e_derivs" not in line:
            this_iter_ederivs.append(float(line.split()[0]))
        elif "</e_derivs" in line:
            collecting_eDerivs = False
            if (len(this_iter_ederivs)>0):
                all_iters_ederivs.append(this_iter_ederivs)
            this_iter_ederivs = []

        if collecting_lDerivs and "</lderivs" not in line:
            this_iter_lderivs.append(float(line.split()[0]))
        elif "</lderivs" in line:
            collecting_lDerivs = False
            if (len(this_iter_lderivs)>0):
                all_iters_lderivs.append(this_iter_lderivs)
            this_iter_lderivs = []
f.close()

# Each iteration is its own plot, and the x axis is the parameters

print("Energy derivatives, normalized")
normalized_ederivs = [norm_param(one_iter_ederivs) for one_iter_ederivs in all_iters_ederivs]
print("Lagrangian derivatives, normalized")
normalized_lderivs = [norm_param(one_iter_lderivs) for one_iter_lderivs in all_iters_lderivs]

num_plots = 3

for i in range(num_plots):
    plt.figure(i)
    plt.plot(np.arange(len(normalized_ederivs[i])), normalized_ederivs[i], label="Energy Derivatives")
    plt.plot(np.arange(len(normalized_lderivs[i])), normalized_lderivs[i], label="Lagrangian Derivatives")
    plt.legend()

plt.show()

