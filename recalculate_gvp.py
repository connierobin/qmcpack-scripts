import os, sys
import matplotlib.pyplot as plt
import numpy as np
import statistics
import math

##################
# SKIM HISTORIES #
##################

def skim_histories():
    print(f"Skimming history vectors")
    # Vectors to hold all samples and iterations, as well as one iteration or one sample
    # that is currently being worked with
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

    collecting_le_der = False
    collecting_der_rat = False
    collecting_lev = False
    collecting_vgs = False

    # Input is the raw QMC Output file
    filename = sys.argv[1]
    print(filename)

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
    if False:
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
    return (all_iters_der_rat, all_iters_le_der, all_iters_lev, all_iters_vgs)


# NOTE: the length of the third dimension for all_iters_le_der and all_iters_der_rat is equal to 
# the number of parameters. The first value is NOT a dud. The first value of der_rat_samp is always
# 1.0 because there is one extra degree of freedom due to normalization. 
# BUT the first CI parameter when doing an excited state is the ground state filling, and der_rat should
# be 1 for this also because that parameter is ZERO! 

######################
# CALCULATE AVERAGES #
######################

# TODO: Do some of these need to be k + 1 instead of k?

def calculate_averages(num_iters, num_samples, num_params, all_iters_der_rat, all_iters_le_der, all_iters_lev, all_iters_vgs):
    # Vectors to hold the averages / sums for all iterations
    avg_le_der_history = []     # Two dimensional: iterations, parameters
    avg_der_rat_history = []    # Two dimensional: iterations, parameters
    avg_lev_history = []        # One dimensional: iterations
    sum_vgs_history = []        # One dimensional: iterations

    # Initialize vectors with zeros
    avg_lev_history = np.zeros(len(all_iters_lev))
    sum_vgs_history = np.zeros(len(all_iters_vgs))

    num_params = len(all_iters_der_rat[0][0])

    # Loop over iterations
    for i in range(num_iters):
        print(f"Calculating averages for iteration {i}")
        # Temporary vectors to hold the averages for all parameters. Once filled, append to main vector.
        avg_le_der_samp = np.zeros(len(all_iters_le_der[i][0]))
        avg_der_rat_samp = np.zeros(len(all_iters_le_der[i][0]))
        # Loop over samples for this iteration
        for j in range(num_samples):
            sum_vgs_history[i] += all_iters_vgs[i][j]
            avg_lev_history[i] += all_iters_lev[i][j] * all_iters_vgs[i][j]     # TODO: is multiplying by vgs correct here??
            # Loop over parameters for this sample for this iteration
            for k in range(num_params):
                avg_le_der_samp[k] += all_iters_le_der[i][j][k] * all_iters_vgs[i][j]
                avg_der_rat_samp[k] += all_iters_der_rat[i][j][k] * all_iters_vgs[i][j]
        # Complete the calculation of taking the average
        avg_lev_history[i] = avg_lev_history[i] / sum_vgs_history[i]
        for k in range(len(all_iters_le_der[i][j])):
            avg_le_der_samp[k] = avg_le_der_samp[k] / sum_vgs_history[i]
            avg_der_rat_samp[k] = avg_der_rat_samp[k] / sum_vgs_history[i]
        avg_le_der_history.append(avg_le_der_samp)
        avg_der_rat_history.append(avg_der_rat_samp)
    if False:
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
    return (avg_der_rat_history, avg_le_der_history, avg_lev_history, sum_vgs_history)


#####################################################
# CALCULATE GVP AND GVP DERIVATIVES  -- QMCPACK WAY #
#####################################################

# avg_le_der_history[i][j] is the i'th iteration at the j'th parameter
# avg_der_rat_history[i][j] is the i'th iteration at the j'th parameter
# avg_lev_history[i] is the i'th iteration of the energy
# sum_vgs_history[i] is the i'th iteration of the some of vgs values

# all_iters_der_rat[i][j][k] is the i'th iteration, the j'th sample, the k'th parameter
# all_iters_le_der[i][j][k] is the i'th iteration, the j'thsample, the k'th parameter
# all_iters_lev[i][j] is the i'th iteration at the j'th sample
# all_iters_vgs[i][j] is the i'th iteration at the j'th sample

# LAGRANGIAN

def calculate_qmcpack_lagrangian(num_iters, num_params, avg_le_der_history, avg_der_rat_history, avg_lev_history):
    all_iters_qmcpack_e_derivs = []
    all_iters_qmcpack_target = []
    this_iter_e_derivs = []
    this_iter_target = 0.0

    for i in range(num_iters):
        print(f"Calculating qmcpack lagrangian for iteration {i}")
        for j in range(num_params):
            this_iter_e_derivs.append(2 * avg_le_der_history[i][j] - 2 * avg_lev_history[i] * avg_der_rat_history[i][j])
            this_iter_target += this_iter_e_derivs[j] ** 2
        all_iters_qmcpack_e_derivs.append(this_iter_e_derivs)
        this_iter_e_derivs = []
        all_iters_qmcpack_target.append(this_iter_target)
        this_iter_target = 0.0
    
    return (all_iters_qmcpack_e_derivs, all_iters_qmcpack_target)

# DERIVATIVE OF LAGRANGIAN

def calculate_qmcpack_lagrangian_derivatives(num_iters, num_samples, num_params, all_iters_qmcpack_e_derivs, avg_der_rat_history, avg_lev_history, sum_vgs_history, all_iters_le_der, all_iters_der_rat, all_iters_vgs):
    all_iters_qmcpack_term1 = []
    all_iters_qmcpack_term2 = []
    all_iters_qmcpack_term3 = []
    all_iters_qmcpack_term4 = []

    all_iters_qmcpack_overlap = []
    all_iters_qmcpack_hamiltonian = []
    all_iters_qmcpack_e_deriv_derivs = []
    all_iters_qmcpack_lderivs = []

    # Terms 3 and 4

    this_iter_term3 = np.zeros(num_params)
    this_iter_term4 = np.zeros(num_params)
    this_iter_qmcpack_e_deriv_derivs = np.zeros((num_params, num_params))
    for i in range(num_iters):
        for l in range(num_params):
            for k in range(num_params):
                this_iter_term3[l] += -4 * all_iters_qmcpack_e_derivs[i][l] * avg_der_rat_history[i][k] * all_iters_qmcpack_e_derivs[i][k]
                this_iter_term4[l] += -4 * all_iters_qmcpack_e_derivs[i][k] * avg_der_rat_history[i][l] * all_iters_qmcpack_e_derivs[i][k]
                this_iter_qmcpack_e_deriv_derivs[l][k] += -2 * all_iters_qmcpack_e_derivs[i][l] * avg_der_rat_history[i][k]
                this_iter_qmcpack_e_deriv_derivs[l][k] += -2 * all_iters_qmcpack_e_derivs[i][k] * avg_der_rat_history[i][l]
        all_iters_qmcpack_term3.append(this_iter_term3)
        this_iter_term3 = np.zeros(num_params)
        all_iters_qmcpack_term4.append(this_iter_term4)
        this_iter_term4 = np.zeros(num_params)
        all_iters_qmcpack_e_deriv_derivs.append(this_iter_qmcpack_e_deriv_derivs)
        this_iter_qmcpack_e_deriv_derivs = np.zeros((num_params, num_params))


    # Terms 1 and 2

    term1_part = 0.0
    term2_part = 0.0

    this_iter_term1 = np.zeros(num_params)
    this_iter_term2 = np.zeros(num_params)
    this_iter_qmcpack_overlap = np.zeros((num_params, num_params))
    this_iter_qmcpack_hamiltonian = np.zeros((num_params, num_params))
    for i in range(num_iters):
        print(f"Calculating qmcpack lagrangian derivatives for iteration {i}")
        this_iter_qmcpack_e_deriv_derivs = all_iters_qmcpack_e_deriv_derivs[i]
        for j in range(num_samples):
            for k in range(num_params):
                # TODO: do I need to change the all_iters call to be for parameter k + 1?
                term1_part += all_iters_le_der[i][j][k] * all_iters_qmcpack_e_derivs[i][k]
                term2_part += all_iters_der_rat[i][j][k] * all_iters_qmcpack_e_derivs[i][k]
            for l in range(num_params):
                # TODO: do I need l + 1???
                this_iter_term1[l] += 4 * all_iters_vgs[i][j] * all_iters_der_rat[i][j][l] * term1_part / sum_vgs_history[i]
                this_iter_term2[l] += -4 * all_iters_vgs[i][j] * all_iters_der_rat[i][j][l] * avg_lev_history[i] * term2_part / sum_vgs_history[i]
                for k in range(num_params):
                    this_iter_qmcpack_e_deriv_derivs[l][k] += 2 * all_iters_vgs[i][j] * all_iters_der_rat[i][j][l] * all_iters_le_der[i][j][k] / sum_vgs_history[i]
                    this_iter_qmcpack_e_deriv_derivs[l][k] += -2 * all_iters_vgs[i][j] * all_iters_der_rat[i][j][l] * avg_lev_history[i] * all_iters_der_rat[i][j][k] / sum_vgs_history[i]
                    this_iter_qmcpack_overlap[l][k] += all_iters_vgs[i][j] * all_iters_der_rat[i][j][l] * all_iters_der_rat[i][j][k] / sum_vgs_history[i]
                    this_iter_qmcpack_hamiltonian[l][k] += all_iters_vgs[i][j] * all_iters_der_rat[i][j][l] * all_iters_le_der[i][j][k] / sum_vgs_history[i]
            term1_part = 0.0
            term2_part = 0.0
        all_iters_qmcpack_term1.append(this_iter_term1)
        this_iter_term1 = np.zeros(num_params)
        all_iters_qmcpack_term2.append(this_iter_term2)
        this_iter_term2 = np.zeros(num_params)
        all_iters_qmcpack_overlap.append(this_iter_qmcpack_overlap)
        this_iter_qmcpack_overlap = np.zeros((num_params, num_params))
        all_iters_qmcpack_hamiltonian.append(this_iter_qmcpack_hamiltonian)
        this_iter_qmcpack_hamiltonian = np.zeros((num_params, num_params))
        all_iters_qmcpack_e_deriv_derivs[i] = this_iter_qmcpack_e_deriv_derivs
        this_iter_qmcpack_e_deriv_derivs = np.zeros((num_params, num_params))

    # Combine Terms 1, 2, 3, and 4

    this_iter_lderivs = np.zeros(num_params)
 
    for i in range(num_iters):
        for k in range(num_params):
            this_iter_lderivs[k] += all_iters_qmcpack_term1[i][k]
            this_iter_lderivs[k] += all_iters_qmcpack_term2[i][k]
            this_iter_lderivs[k] += all_iters_qmcpack_term3[i][k]
            this_iter_lderivs[k] += all_iters_qmcpack_term4[i][k]
        all_iters_qmcpack_lderivs.append(this_iter_lderivs)
        this_iter_lderivs = np.zeros(num_params)

    return ((all_iters_qmcpack_term1, all_iters_qmcpack_term2, all_iters_qmcpack_term3, all_iters_qmcpack_term4), all_iters_qmcpack_lderivs, all_iters_qmcpack_overlap, all_iters_qmcpack_hamiltonian, all_iters_qmcpack_e_deriv_derivs)

##################################################
# CALCULATE GVP AND DERIVATIVES -- ALGEBRAIC WAY #
##################################################

def calculate_algebraic_ederivs(num_iters, num_params, avg_der_rat_history, avg_le_der_history, avg_lev_history):
    all_iters_algebraic_e_derivs = []
    # loop through all iterations
    for i in range(num_iters):
        print(f"Calculating algebraic energy derivatives for iteration {i}")
        E = avg_lev_history[i]
        B = avg_le_der_history[i]
        D = avg_der_rat_history[i]
        dEdck = np.zeros(num_params)    # energy derivatives for all parameters
        # loop through all parameters
        for k in range(num_params):
            # formula for energy derivatives
            dEdck[k] = 2 * B[k] - 2 * E * D[k]
        all_iters_algebraic_e_derivs.append(dEdck)
    return all_iters_algebraic_e_derivs

def calculate_algebraic_lagrangian(num_iters, num_params, all_iters_algebraic_e_derivs):
    all_iters_algebraic_target = []
    # loop through all iterations
    for i in range(num_iters):
        print(f"Calculating algebraic lagrangian for iteration {i}")
        dEdck = all_iters_algebraic_e_derivs[i]
        target = 0.0
        # loop through all parameters to sum up contributions to target function
        for k in range(num_params):
            target += dEdck[k] * dEdck[k]
        all_iters_algebraic_target.append(target)
    return all_iters_algebraic_target

def calculate_algebraic_lagrangian_derivatives(num_iters, num_params, all_iters_hamiltonian, all_iters_overlap, all_iters_algebraic_e_derivs, avg_der_rat_history, avg_lev_history):
    all_iters_algebraic_e_deriv_derivs = []
    all_iters_algebraic_lderivs = []

    # Calculate energy second derivatives
    for i in range(num_iters):
        print(f"Calculating algebraic energy second derivatives for iteration {i}")
        this_iter_e_deriv_derivs = np.zeros((num_params, num_params))
        hamiltonian = all_iters_hamiltonian[i]
        overlap = all_iters_overlap[i]
        avgD = avg_der_rat_history[i]
        E = avg_lev_history[i]
        dEdck = all_iters_algebraic_e_derivs[i]
        for k in range(num_params):
            for l in range(num_params):
                this_iter_e_deriv_derivs[l][k] = 2 * hamiltonian[l][k] - 2 * E * overlap[l][k] - 2 * dEdck[l] * avgD[k] - 2 * avgD[l] * dEdck[k]
        all_iters_algebraic_e_deriv_derivs.append(this_iter_e_deriv_derivs)

    # Calculate lagrangian derivatives
    for i in range(num_iters):
        print(f"Calculating algebraic lagrangian derivatives for iteration {i}")
        this_iter_lderivs = np.zeros(num_params)
        e_deriv_derivs = all_iters_algebraic_e_deriv_derivs[i]
        e_derivs = all_iters_algebraic_e_derivs[i]
        for k in range(num_params):
            for l in range(num_params):
                this_iter_lderivs[l] += 2 * e_deriv_derivs[l][k] * e_derivs[k]
        all_iters_algebraic_lderivs.append(this_iter_lderivs)
    
    return (all_iters_algebraic_e_deriv_derivs, all_iters_algebraic_lderivs)

def calculate_algebraic_lagrangian_derivative_terms(num_iters, num_params, all_iters_hamiltonian, all_iters_overlap, all_iters_algebraic_e_derivs, avg_lev_history, avg_der_rat_history):
    all_iters_algebraic_term1 = []
    all_iters_algebraic_term2 = []
    all_iters_algebraic_term3 = []
    all_iters_algebraic_term4 = []
    all_iters_algebraic_terms_sum = []

    for i in range(num_iters):
        this_iter_algebraic_term1 = np.zeros(num_params)
        this_iter_algebraic_term2 = np.zeros(num_params)
        this_iter_algebraic_term3 = np.zeros(num_params)
        this_iter_algebraic_term4 = np.zeros(num_params)
        for k in range(num_params):
            for l in range(num_params):
                this_iter_algebraic_term1[l] += 4 * all_iters_hamiltonian[i][l][k] * all_iters_algebraic_e_derivs[i][k]
                this_iter_algebraic_term2[l] += -4 * avg_lev_history[i] * all_iters_overlap[i][l][k] * all_iters_algebraic_e_derivs[i][k]
                this_iter_algebraic_term3[l] += -4 * all_iters_algebraic_e_derivs[i][l] * avg_der_rat_history[i][k] * all_iters_algebraic_e_derivs[i][k]
                this_iter_algebraic_term4[l] += -4 * avg_der_rat_history[i][l] * all_iters_algebraic_e_derivs[i][k] * all_iters_algebraic_e_derivs[i][k]
        all_iters_algebraic_term1.append(this_iter_algebraic_term1)
        all_iters_algebraic_term2.append(this_iter_algebraic_term2)
        all_iters_algebraic_term3.append(this_iter_algebraic_term3)
        all_iters_algebraic_term4.append(this_iter_algebraic_term4)

    for i in range(num_iters):
        this_iter_algebraic_terms_sum = np.zeros(num_params)
        for k in range(num_params):
            this_iter_algebraic_terms_sum[k] = all_iters_algebraic_term1[i][k] + all_iters_algebraic_term2[i][k] + all_iters_algebraic_term3[i][k] + all_iters_algebraic_term4[i][k]
        all_iters_algebraic_terms_sum.append(this_iter_algebraic_terms_sum)

    return ((all_iters_algebraic_term1, all_iters_algebraic_term2, all_iters_algebraic_term3, all_iters_algebraic_term4), all_iters_algebraic_terms_sum)

##############################################
# CALCULATE OVERLAP AND HAMILTONIAN MATRICES #
##############################################

def calculate_hamiltonian_and_overlap(num_iters, num_samples, num_params, all_iters_der_rat, all_iters_le_der, all_iters_vgs, sum_vgs_history):
    all_iters_overlap = []
    all_iters_hamiltonian = []

    # Calculate for each iteration
    for i in range(num_iters):
        print(f"Calculating overlap and hamiltonian for iteration {i}")
        # Reset the matrices for this iteration so we can sum up starting from zero
        this_iter_overlap = np.zeros((num_params, num_params))
        this_iter_hamiltonian = np.zeros((num_params, num_params))
        vgs_sum = sum_vgs_history[i]    # normalization factor
        # Add up the contribution from each sample
        for j in range(num_samples):
            D = all_iters_der_rat[i][j] # phi_k / psi for iteration i sample j
            B = all_iters_le_der[i][j]  # H phi_k / psi for iteration i sample j
            vgs = all_iters_vgs[i][j]   # value guiding squared for iteration i sample j
            # iterate through the matrix elements
            for k in range(num_params):
                for l in range(num_params):
                    this_iter_overlap[k][l] += vgs * D[k] * D[l]
                    this_iter_hamiltonian[k][l] += vgs * D[k] * B[l]
        # Normalize, element by element
        for k in range(num_params):
            for l in range(num_params):
                this_iter_overlap[k][l] = this_iter_overlap[k][l] / vgs_sum
                this_iter_hamiltonian[k][l] = this_iter_hamiltonian[k][l] / vgs_sum
        # Store the results
        all_iters_overlap.append(this_iter_overlap)
        all_iters_hamiltonian.append(this_iter_hamiltonian)

    return (all_iters_overlap, all_iters_hamiltonian)

#################################################
# COMPARE QMCPACK NUMBERS AND ALGEBRAIC NUMBERS #
#################################################

def compare_e_derivs(num_iters, num_params, all_iters_qmcpack_e_derivs, all_iters_algebraic_e_derivs):
    print("Comparing energy derivatives...")
    for i in range(num_iters):
        print(f"Iteration {i}")
        num_mismatch = 0    # tally up the number of elements that don't match
        avg_mismatch = 0.0   # average distance between the numbers
        for k in range(num_params):
            if not math.isclose(all_iters_qmcpack_e_derivs[i][k], all_iters_algebraic_e_derivs[i][k]):
                num_mismatch += 1
                avg_mismatch += np.abs(all_iters_qmcpack_e_derivs[i][k] - all_iters_algebraic_e_derivs[i][k])
        if num_mismatch > 0:
            avg_mismatch = avg_mismatch / num_mismatch
            print(f"Number of mismatches: {num_mismatch}")
            print(f"Average distance between values: {avg_mismatch}")
        else:
            print(f"No mismatch found!")

def compare_target(num_iters, all_iters_qmcpack_target, all_iters_algebraic_target):
    print("Comparing Lagrangians...")
    for i in range(num_iters):
        print(f"Iteration {i}")
        if not math.isclose(all_iters_qmcpack_target[i], all_iters_algebraic_target[i]):
            print(f"Mismatch found! Distance between values: {np.abs(all_iters_qmcpack_target[i] - all_iters_algebraic_target[i])}")
        else:
            print(f"No mismatch found!")

def compare_lderivs(num_iters, num_params, all_iters_qmcpack_lderivs, all_iters_algebraic_lderivs, all_iters_algebraic_terms_sum):
    print("Comparing Lagrangian derivatives between qmcpack and algebraic...")
    for i in range(num_iters):
        print(f"Iteration {i}")
        num_mismatch = 0    # tally up the number of elements that don't match
        avg_mismatch = 0.0   # average distance between the numbers
        for k in range(num_params):
            if not math.isclose(all_iters_qmcpack_lderivs[i][k], all_iters_algebraic_lderivs[i][k]):
                num_mismatch += 1
                avg_mismatch += np.abs(all_iters_qmcpack_lderivs[i][k] - all_iters_algebraic_lderivs[i][k])
                # print(f"Mismatch found! QMCPACK number: {all_iters_qmcpack_lderivs[i][k]}   algebraic number: {all_iters_algebraic_lderivs[i][k]}")
        if num_mismatch > 0:
            avg_mismatch = avg_mismatch / num_mismatch
            print(f"Number of mismatches: {num_mismatch}")
            print(f"Average distance between values: {avg_mismatch}")
        else:
            print(f"No mismatch found!")

    print("Comparing Lagrangian derivatives between algebraic and algebraic terms sum...")
    for i in range(num_iters):
        print(f"Iteration {i}")
        num_mismatch = 0    # tally up the number of elements that don't match
        avg_mismatch = 0.0   # average distance between the numbers
        for k in range(num_params):
            if not math.isclose(all_iters_algebraic_lderivs[i][k], all_iters_algebraic_terms_sum[i][k]):
                num_mismatch += 1
                avg_mismatch += np.abs(all_iters_qmcpack_lderivs[i][k] - all_iters_algebraic_lderivs[i][k])
                # print(f"Mismatch found! QMCPACK number: {all_iters_qmcpack_lderivs[i][k]}   algebraic number: {all_iters_algebraic_lderivs[i][k]}")
        if num_mismatch > 0:
            avg_mismatch = avg_mismatch / num_mismatch
            print(f"Number of mismatches: {num_mismatch}")
            print(f"Average distance between values: {avg_mismatch}")
        else:
            print(f"No mismatch found!")

def compare_lderiv_terms(num_iters, num_params, all_iters_qmcpack_lderiv_terms, all_iters_algebraic_lderiv_terms):
    print("Comparing Lagrangian derivative terms...")
    
    (all_iters_algebraic_term1, all_iters_algebraic_term2, all_iters_algebraic_term3, all_iters_algebraic_term4) = all_iters_algebraic_lderiv_terms
    (all_iters_qmcpack_term1, all_iters_qmcpack_term2, all_iters_qmcpack_term3, all_iters_qmcpack_term4) = all_iters_qmcpack_lderiv_terms

    for i in range(num_iters):
        print(f"Iteration {i}: Comparing Term 1...")
        compare_lderiv_term(num_params, all_iters_qmcpack_term1[i], all_iters_algebraic_term1[i])
        print(f"Iteration {i}: Comparing Term 2...")
        compare_lderiv_term(num_params, all_iters_qmcpack_term2[i], all_iters_algebraic_term2[i])
        print(f"Iteration {i}: Comparing Term 3...")
        compare_lderiv_term(num_params, all_iters_qmcpack_term3[i], all_iters_algebraic_term3[i])
        print(f"Iteration {i}: Comparing Term 4...")
        compare_lderiv_term(num_params, all_iters_qmcpack_term4[i], all_iters_algebraic_term4[i])

def compare_lderiv_term(num_params, this_iter_qmcpack_lderiv_term, this_iter_algebraic_lderiv_term):    
    num_mismatch = 0    # tally up the number of elements that don't match
    avg_mismatch = 0.0   # average distance between the numbers
    for k in range(num_params):
        if not math.isclose(this_iter_qmcpack_lderiv_term[k], this_iter_algebraic_lderiv_term[k]):
            num_mismatch += 1
            avg_mismatch += np.abs(this_iter_qmcpack_lderiv_term[k] - this_iter_algebraic_lderiv_term[k])
            print(f"Mismatch found! QMCPACK number: {this_iter_qmcpack_lderiv_term[k]}   algebraic number: {this_iter_algebraic_lderiv_term[k]}")
    if num_mismatch > 0:
        avg_mismatch = avg_mismatch / num_mismatch
        print(f"Number of mismatches: {num_mismatch}")
        print(f"Average distance between values: {avg_mismatch}")
    else:
        print(f"No mismatch found!")

def compare_e_deriv_derivs(num_iters, num_params, all_iters_qmcpack_e_deriv_derivs, all_iters_algebraic_e_deriv_derivs):
    print("Comparing energy second derivatives...")
    for i in range(num_iters):
        print(f"Iteration {i}")
        num_mismatch = 0    # tally up the number of elements that don't match
        avg_mismatch = 0.0   # average distance between the numbers
        for k in range(num_params):
            for l in range(num_params):
                if not math.isclose(all_iters_qmcpack_e_deriv_derivs[i][l][k], all_iters_algebraic_e_deriv_derivs[i][l][k]):
                    num_mismatch += 1
                    avg_mismatch += np.abs(all_iters_qmcpack_e_deriv_derivs[i][l][k] - all_iters_algebraic_e_deriv_derivs[i][l][k])
        if num_mismatch > 0:
            avg_mismatch = avg_mismatch / num_mismatch
            print(f"Number of mismatches: {num_mismatch}")
            print(f"Average distance between values: {avg_mismatch}")
        else:
            print(f"No mismatch found!")

def compare_hamiltonian_and_overlap(num_iters, num_params, all_iters_overlap, all_iters_qmcpack_overlap, all_iters_hamiltonian, all_iters_qmcpack_hamiltonian):
    print("Comparing overlap matrix...")
    for i in range(num_iters):
        print(f"Iteration {i}")
        num_mismatch = 0    # tally up the number of elements that don't match
        avg_mismatch = 0.0   # average distance between the numbers
        for k in range(num_params):
            for l in range(num_params):
                if not math.isclose(all_iters_qmcpack_overlap[i][l][k], all_iters_overlap[i][l][k]):
                    num_mismatch += 1
                    avg_mismatch += np.abs(all_iters_qmcpack_overlap[i][l][k] - all_iters_overlap[i][l][k])
        if num_mismatch > 0:
            avg_mismatch = avg_mismatch / num_mismatch
            print(f"Number of mismatches: {num_mismatch}")
            print(f"Average distance between values: {avg_mismatch}")
        else:
            print(f"No mismatch found!")

    print("Comparing hamiltonian matrix...")
    for i in range(num_iters):
        print(f"Iteration {i}")
        num_mismatch = 0    # tally up the number of elements that don't match
        avg_mismatch = 0.0   # average distance between the numbers
        for k in range(num_params):
            for l in range(num_params):
                if not math.isclose(all_iters_qmcpack_hamiltonian[i][l][k], all_iters_hamiltonian[i][l][k]):
                    num_mismatch += 1
                    avg_mismatch += np.abs(all_iters_qmcpack_hamiltonian[i][l][k] - all_iters_hamiltonian[i][l][k])
        if num_mismatch > 0:
            avg_mismatch = avg_mismatch / num_mismatch
            print(f"Number of mismatches: {num_mismatch}")
            print(f"Average distance between values: {avg_mismatch}")
        else:
            print(f"No mismatch found!")

def compare_qmcpack_algebraic(num_iters, num_params, all_iters_qmcpack_e_derivs, all_iters_algebraic_e_derivs, all_iters_qmcpack_target, all_iters_algebraic_target, all_iters_qmcpack_lderiv_terms, all_iters_algebraic_lderiv_terms, all_iters_qmcpack_lderivs, all_iters_algebraic_lderivs, all_iters_algebraic_terms_sum, all_iters_qmcpack_e_deriv_derivs, all_iters_algebraic_e_deriv_derivs, all_iters_overlap, all_iters_qmcpack_overlap, all_iters_hamiltonian, all_iters_qmcpack_hamiltonian):
    # TODO: make sure that num_iters and num_params is consistent across all results. use array sizes to compare
    # ALL results should be normalized to have a true comparison. 
    compare_e_derivs(num_iters, num_params, all_iters_qmcpack_e_derivs, all_iters_algebraic_e_derivs)
    compare_target(num_iters, all_iters_qmcpack_target, all_iters_algebraic_target)
    compare_lderivs(num_iters, num_params, all_iters_qmcpack_lderivs, all_iters_algebraic_lderivs, all_iters_algebraic_terms_sum)
    compare_lderiv_terms(num_iters, num_params, all_iters_qmcpack_lderiv_terms, all_iters_algebraic_lderiv_terms)
    compare_e_deriv_derivs(num_iters, num_params, all_iters_qmcpack_e_deriv_derivs, all_iters_algebraic_e_deriv_derivs)
    compare_hamiltonian_and_overlap(num_iters, num_params, all_iters_overlap, all_iters_qmcpack_overlap, all_iters_hamiltonian, all_iters_qmcpack_hamiltonian)

#####################
# PRINTER FUNCTIONS #
#####################

def print_all_elements(num_iters, num_samples, num_params, all_iters_der_rat, all_iters_le_der, all_iters_lev, all_iters_vgs):
    for i in range(num_iters):
        print(f"\nIteration {i}")
        print(f"\nder_rat_samp")
        for j in range(num_samples):
            print(f"Sample {j}")
            this_row = []
            for k in range(num_params):
                this_row.append(all_iters_der_rat[i][j][k])
            print("\t".join([str(elem) for elem in this_row]))
        print(f"\nle_der_samp")
        for j in range(num_samples):
            print(f"Sample {j}")
            this_row = []
            for k in range(num_params):
                this_row.append(all_iters_le_der[i][j][k])
            print("\t".join([str(elem) for elem in this_row]))
        print(f"\nvgs")
        for j in range(num_samples):
            print(f"Sample {j}")
            print(all_iters_vgs[i][j])
        print(f"\nEnergy")
        for j in range(num_samples):
            print(f"Sample {j}")
            print(all_iters_lev[i][j])

def print_overlap_and_hamiltonian(num_iters, num_params, all_iters_overlap, all_iters_qmcpack_overlap, all_iters_hamiltonian, all_iters_qmcpack_hamiltonian):
    print("\nAlgebraic Overlap: ")
    for i in range(num_iters):
        for j in range(num_params):
            this_row = []
            for k in range(num_params):
                this_row.append(all_iters_overlap[i][j][k])
            print("\t".join([str(elem) for elem in this_row]))
    print("end")
    print("\nQMCPACK Overlap: ")
    for i in range(num_iters):
        for j in range(num_params):
            this_row = []
            for k in range(num_params):
                this_row.append(all_iters_qmcpack_overlap[i][j][k])
            print("\t".join([str(elem) for elem in this_row]))
    print("end")
    print("\nAlgebraic Hamiltonian: ")
    for i in range(num_iters):
        for j in range(num_params):
            this_row = []
            for k in range(num_params):
                this_row.append(all_iters_hamiltonian[i][j][k])
            print("\t".join([str(elem) for elem in this_row]))
    print("end")
    print("\nQMCPACK Hamiltonian: ")
    for i in range(num_iters):
        for j in range(num_params):
            this_row = []
            for k in range(num_params):
                this_row.append(all_iters_qmcpack_hamiltonian[i][j][k])
            print("\t".join([str(elem) for elem in this_row]))
    print("end")

def print_lderivs(num_iters, num_params, all_iters_algebraic_lderivs, all_iters_qmcpack_lderivs):
    print("\nAlgebraic Lagrangian derivatives: ")
    lderiv_magnitude = 0.0
    for i in range(num_iters):
        for j in range(num_params):
            lderiv_magnitude += all_iters_algebraic_lderivs[i][j] ** 2
            print(f"{all_iters_algebraic_lderivs[i][j]}\n")
    print("end")
    lderiv_magnitude = np.sqrt(lderiv_magnitude)
    print(f"\nAlgebraic Lagrangian derivative magnitude: {lderiv_magnitude}")

    print("\QMCPACK Lagrangian derivatives: ")
    lderiv_magnitude = 0.0
    for i in range(num_iters):
        for j in range(num_params):
            lderiv_magnitude += all_iters_qmcpack_lderivs[i][j] ** 2
            print(f"{all_iters_qmcpack_lderivs[i][j]}\n")
    print("end")
    lderiv_magnitude = np.sqrt(lderiv_magnitude)
    print(f"\QMCPACK Lagrangian derivative magnitude: {lderiv_magnitude}")

########
# MAIN #
########

def main():
    (all_iters_der_rat, all_iters_le_der, all_iters_lev, all_iters_vgs)                                                             = skim_histories()
    num_samples = len(all_iters_der_rat[0])     # CHANGE this line to be a small number to make it take less time
    num_iters = len(all_iters_der_rat)       # CHANGE this line to be a small number to make it take less time
    num_params = len(all_iters_der_rat[0][0])
    (avg_der_rat_history, avg_le_der_history, avg_lev_history, sum_vgs_history)                                                     = calculate_averages(num_iters, num_samples, num_params, all_iters_der_rat, all_iters_le_der, all_iters_lev, all_iters_vgs)
    (all_iters_qmcpack_e_derivs, all_iters_qmcpack_target)                                                                          = calculate_qmcpack_lagrangian(num_iters, num_params, avg_le_der_history, avg_der_rat_history, avg_lev_history)
    (all_iters_qmcpack_lderiv_terms, all_iters_qmcpack_lderivs, all_iters_qmcpack_overlap, all_iters_qmcpack_hamiltonian, all_iters_qmcpack_e_deriv_derivs)    = calculate_qmcpack_lagrangian_derivatives(num_iters, num_samples, num_params, all_iters_qmcpack_e_derivs, avg_der_rat_history, avg_lev_history, sum_vgs_history, all_iters_le_der, all_iters_der_rat, all_iters_vgs)
    (all_iters_overlap, all_iters_hamiltonian)                                                                                      = calculate_hamiltonian_and_overlap(num_iters, num_samples, num_params, all_iters_der_rat, all_iters_le_der, all_iters_vgs, sum_vgs_history)
    all_iters_algebraic_e_derivs                                                                                                    = calculate_algebraic_ederivs(num_iters, num_params, avg_der_rat_history, avg_le_der_history, avg_lev_history)
    all_iters_algebraic_target                                                                                                      = calculate_algebraic_lagrangian(num_iters, num_params, all_iters_algebraic_e_derivs)
    (all_iters_algebraic_e_deriv_derivs, all_iters_algebraic_lderivs)                                                               = calculate_algebraic_lagrangian_derivatives(num_iters, num_params, all_iters_hamiltonian, all_iters_overlap, all_iters_algebraic_e_derivs, avg_der_rat_history, avg_lev_history)
    (all_iters_algebraic_lderiv_terms, all_iters_algebraic_terms_sum)                                                               = calculate_algebraic_lagrangian_derivative_terms(num_iters, num_params, all_iters_hamiltonian, all_iters_overlap, all_iters_algebraic_e_derivs, avg_lev_history, avg_der_rat_history)
    compare_qmcpack_algebraic(num_iters, num_params, all_iters_qmcpack_e_derivs, all_iters_algebraic_e_derivs, all_iters_qmcpack_target, all_iters_algebraic_target, all_iters_qmcpack_lderiv_terms, all_iters_algebraic_lderiv_terms, all_iters_qmcpack_lderivs, all_iters_algebraic_lderivs, all_iters_algebraic_terms_sum, all_iters_qmcpack_e_deriv_derivs, all_iters_algebraic_e_deriv_derivs, all_iters_overlap, all_iters_qmcpack_overlap, all_iters_hamiltonian, all_iters_qmcpack_hamiltonian)
    # print_all_elements(num_iters, num_samples, num_params, all_iters_der_rat, all_iters_le_der, all_iters_lev, all_iters_vgs)
    print_lderivs(num_iters, num_params, all_iters_algebraic_lderivs, all_iters_qmcpack_lderivs)
    print_overlap_and_hamiltonian(num_iters, num_params, all_iters_overlap, all_iters_qmcpack_overlap, all_iters_hamiltonian, all_iters_qmcpack_hamiltonian)
    # TODO: make these functions to see if my calculations are the same as what is happening within qmcpack
    # skim_gvp_terms()
    # compare_skim_gvp_to_calculated()
    # TODO: correct both the algebraic and qmcpack versions to correctly exclude the first number to make sure that can still line up

if __name__ == "__main__":
    main()