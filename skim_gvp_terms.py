import os, sys
import matplotlib.pyplot as plt
import numpy as np
import statistics

#TODO:
# Graph the results non-normalized

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

collecting_term1 = False
collecting_term2 = False
collecting_term3 = False
collecting_term4 = False

this_iter_term1 = []
all_iters_term1=[]
this_iter_term2 = []
all_iters_term2=[]
this_iter_term3 = []
all_iters_term3=[]
this_iter_term4 = []
all_iters_term4=[]

f=open(filename,'r')
for line in f:
    # Collecting set to true upon finding an <eDerivs> tag
    # Stop collecting when we hit the end of the variable section
    if not collecting_term1 and not collecting_term2 and not collecting_term3 and not collecting_term4:
        if "<term_1" in line:
            collecting_term1 = True
        if "<term_2" in line:
            collecting_term2 = True
        if "<term_3" in line:
            collecting_term3 = True
        if "<term_4" in line:
            collecting_term4 = True
    else:
        if collecting_term1 and "</term_1" not in line:
            this_iter_term1.append(float(line.split()[0]))
        elif "</term_1" in line:
            collecting_term1 = False
            if (len(this_iter_term1)>0):
                all_iters_term1.append(this_iter_term1)
            this_iter_term1 = []

        if collecting_term2 and "</term_2" not in line:
            this_iter_term2.append(float(line.split()[0]))
        elif "</term_2" in line:
            collecting_term2 = False
            if (len(this_iter_term2)>0):
                all_iters_term2.append(this_iter_term2)
            this_iter_term2 = []

        if collecting_term3 and "</term_3" not in line:
            this_iter_term3.append(float(line.split()[0]))
        elif "</term_3" in line:
            collecting_term3 = False
            if (len(this_iter_term3)>0):
                all_iters_term3.append(this_iter_term3)
            this_iter_term3 = []

        if collecting_term4 and "</term_4" not in line:
            this_iter_term4.append(float(line.split()[0]))
        elif "</term_4" in line:
            collecting_term4 = False
            if (len(this_iter_term4)>0):
                all_iters_term4.append(this_iter_term4)
            this_iter_term4 = []
f.close()

# Each iteration is its own plot, and the x axis is the parameters

print("Term 1 values, normalized")
normalized_term1 = [norm_param(one_iter_term) for one_iter_term in all_iters_term1]
print("Term 2 values, normalized")
normalized_term2 = [norm_param(one_iter_term) for one_iter_term in all_iters_term2]
print("Term 3 values, normalized")
normalized_term3 = [norm_param(one_iter_term) for one_iter_term in all_iters_term3]
print("Term 4 values, normalized")
normalized_term4 = [norm_param(one_iter_term) for one_iter_term in all_iters_term4]

num_plots = 3

for i in range(num_plots):
    plt.figure(i)
    plt.plot(np.arange(len(normalized_term1[i])), normalized_term1[i], label="Term 1")
    plt.plot(np.arange(len(normalized_term2[i])), normalized_term2[i], label="Term 2")
    plt.plot(np.arange(len(normalized_term3[i])), normalized_term3[i], label="Term 3")
    plt.plot(np.arange(len(normalized_term4[i])), normalized_term4[i], label="Term 4")
    plt.legend()

plt.show()

