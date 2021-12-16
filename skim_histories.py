import os, sys
import matplotlib.pyplot as plt
import numpy as np
import statistics

#TODO:
#skim the history vectors
#store the results
#sotre them for each iteration
    
# Input is the raw QMC Output file
filename = sys.argv[1]
print(filename)

collecting_le_der = False
collecting_der_rat = False
collecting_lev = False
collecting_vgs = False

this_samp_le_der = []
this_iter_le_der = []
all_iters_le_der=[]

this_samp_der_rat = []
this_iter_der_rat = []
all_iters_der_rat=[]

this_iter_lev = []
all_iters_lev=[]

this_iter_vgs = []
all_iters_vgs=[]

f=open(filename,'r')
for line in f:
    # Collecting set to true upon finding an <eDerivs> tag
    # Stop collecting when we hit the end of the variable section
    if not collecting_le_der and not collecting_der_rat and not collecting_lev and not collecting_vgs:
        if "<le_der_history" in line:
            collecting_le_der = True
        if "<der_rat_history" in line:
            collecting_der_rat = True
        if "<lev_history" in line:
            collecting_lev = True
        if "<vg_history" in line:
            collecting_vgs = True
    else:
        if collecting_le_der and "</le_der_history" not in line:
            if "</" in line:
                this_iter_le_der.append(this_samp_le_der)
                this_samp_le_der = []
            elif "<" in line:
                pass
            else:
                this_samp_le_der.append(float(line.split()[0]))
        elif "</le_der_history" in line:
            collecting_le_der = False
            if (len(this_iter_le_der)>0):
                all_iters_le_der.append(this_iter_le_der)
            this_iter_le_der = []

        if collecting_der_rat and "</der_rat_history" not in line:
            if "</" in line:
                this_iter_der_rat.append(this_samp_der_rat)
                this_samp_der_rat = []
            elif "<" in line:
                pass
            else:
                this_samp_der_rat.append(float(line.split()[0]))
        elif "</der_rat_history" in line:
            collecting_der_rat = False
            if (len(this_iter_der_rat)>0):
                all_iters_der_rat.append(this_iter_der_rat)
            this_iter_der_rat = []

        if collecting_lev and "</lev_history" not in line:
            this_iter_lev.append(float(line.split()[0]))
        elif "</lev_history" in line:
            collecting_lev = False
            if (len(this_iter_lev)>0):
                all_iters_lev.append(this_iter_lev)
            this_iter_lev = []

        if collecting_vgs and "</vg_history" not in line:
            this_iter_vgs.append(float(line.split()[0]))
        elif "</vg_history" in line:
            collecting_vgs = False
            if (len(this_iter_vgs)>0):
                all_iters_vgs.append(this_iter_vgs)
            this_iter_vgs = []

f.close()

print("\nle_der_history vector dimensions")
print("Number of iterations: ", len(all_iters_le_der))
print("Number of samples: ", len(all_iters_le_der[0]))
print("Number of parameters: ", len(all_iters_le_der[0][0]))
print("\nder_rat_history vector dimensions")
print("Number of iterations: ", len(all_iters_der_rat))
print("Number of samples: ", len(all_iters_der_rat[0]))
print("Number of parameters: ", len(all_iters_der_rat[0][0]))
print("\nlev_history vector dimensions")
print("Number of iterations: ", len(all_iters_lev))
print("Number of samples: ", len(all_iters_lev[0]))
print("\nvgs_history vector dimensions")
print("Number of iterations: ", len(all_iters_vgs))
print("Number of samples: ", len(all_iters_vgs[0]))


# NOTE: the length of the third dimension for all_iters_le_der and all_iters_der_rat is equal to 
# the number of parameters. The first value is NOT a dud. The first value of der_rat_samp is always
# 1.0 because there is one extra degree of freedom due to normalization. 

# CALCULATE AVERAGES

# Vectors to hold the averages / sums for all iterations
avg_le_der_history = []     # Two dimensional: iterations, parameters
avg_der_rat_history = []    # Two dimensional: iterations, parameters
avg_lev_history = []        # One dimensional: iterations
sum_vgs_history = []        # One dimensional: iterations
# Initialize vectors with zeros
avg_lev_history = np.zeros(len(all_iters_lev))
sum_vgs_history = np.zeros(len(all_iters_vgs))

# Loop over iterations
for i in range(len(all_iters_le_der)):
    print(f"Calculating averages for iteration {i}")
    # Temporary vectors to hold the averages for all parameters. Once filled, append to main vector.
    avg_le_der_samp = np.zeros(len(all_iters_le_der[i][0]))
    avg_der_rat_samp = np.zeros(len(all_iters_le_der[i][0]))
    # Loop over samples for this iteration
    for j in range(len(all_iters_le_der[i])):
        sum_vgs_history[i] += all_iters_vgs[i][j]
        avg_lev_history[i] += all_iters_lev[i][j]
        # Loop over parameters for this sample for this iteration
        for k in range(len(all_iters_le_der[i][j])):
            avg_le_der_samp[k] += all_iters_le_der[i][j][k]
            avg_der_rat_samp[k] += all_iters_der_rat[i][j][k]
    # Complete the calculation of taking the average
    avg_lev_history[i] = avg_lev_history[i] / sum_vgs_history[i]
    for k in range(len(all_iters_le_der[i][j])):
        avg_le_der_samp[k] = avg_le_der_samp[k] / sum_vgs_history[i]
        avg_der_rat_samp[k] = avg_der_rat_samp[k] / sum_vgs_history[i]
    avg_le_der_history.append(avg_le_der_samp)
    avg_der_rat_history.append(avg_der_rat_samp)

# Print dimensions to see how we're doing
print("avg_le_der_history dimensions: ")
print("Number of iterations: ", len(avg_le_der_history))
print("Number of parameters: ", len(avg_le_der_history[0]))
print("avg_der_rat_history dimensions: ")
print("Number of iterations: ", len(avg_der_rat_history))
print("Number of parameters: ", len(avg_der_rat_history[0]))
print("avg_lev_history dimensions: ")
print("Number of iterations: ", len(avg_lev_history))
print("sum_vgs_history dimensions: ")
print("Number of iterations: ", len(sum_vgs_history))




    # for (int i = 0; i < num_optimizables; i++)
    # {
    #   avg_le_der_samp_.at(i) = avg_le_der_samp_.at(i) / vgs_sum_;
    #   avg_der_rat_samp_.at(i) = avg_der_rat_samp_.at(i) / vgs_sum_;

    #   e_derivs_.push_back(2 * avg_le_der_samp_.at(i) - 2 * e_avg_ * avg_der_rat_samp_.at(i));
    #   //TODO:Connie later this will need to have the mpi_unbiased_ratio applied to it most likely.
    #   target_avg_ += e_derivs_.at(i) * e_derivs_.at(i);
    # }