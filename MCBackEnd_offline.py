from MCBackend import MCBackend
import json
import os

"""This code simulates lightning strikes and calculates the TDOA for the given network of stations without the user interface"""

# DEFINITIONS 
num_ls = 10000                              # number of lightning strikes to simulate
region_of_network = 'SampleNetwork-USA'     # region of the network to be used
mode = False                                # True: all stations will detect all lightning. It will ignore the models in the pkl files.
                                            # False: considers the dynamics of detectors and the model (pkl) of the station response is used. 

# INPUT FILES
input_stations_file = 'STA_network.json'                        # station network file
curr_dist_parameters_file = 'curr_distribution_parameters.json' # file with the parameters of the current distribution
svm_parameters_negLS = 'svm_peak_neg_300km.pkl'                 # file with the parameters of the SVM model
svm_parameters_posLS = 'svm_peak_pos_300km.pkl'                 # file with the parameters of the SVM model
tdoa_solve_file = 'tdoa_solve.py'                               # file with the TDOA solver

# OUTPUT FILES
output_file_sim = 'LSsimulated_offline.json'                    # output file for simulated lightning strikes
output_file_res = 'LSresults_offline.json'                      # output file for results

# Check if the input files exist
input_files = [input_stations_file, curr_dist_parameters_file, svm_parameters_negLS, svm_parameters_posLS, tdoa_solve_file]
missing_files = [file for file in input_files if not os.path.isfile(file)]

if input_stations_file in missing_files or curr_dist_parameters_file in missing_files:
    print("The following input files do not exist, and it is impossible to continue:")
    for file in missing_files:
        print(file)
else:
    if missing_files:
        print("The following input files do not exist, but the code can continue:")
        for file in missing_files:
            print(file)
    else:
        print("All input files exist.")

    bckend = MCBackend(input_stations_file, curr_dist_parameters_file, svm_parameters_negLS, svm_parameters_posLS)

    with open(input_stations_file, 'r') as json_file:               # load station network
        STA_network = json.load(json_file)
    input_stations = bckend.filter_stations_by_region(region_of_network)        # selects desired network among stations dictionary

    bckend.mode = mode                                              # if False, considers the actual dynamics of detectors
    bckend.add_input_stations(input_stations)                       # add stations to backend for solution
    lstrikes = bckend.commander(num_ls)                             # generate num_ls lightning strikes
    bckend.save(lstrikes,filename=output_file_sim)                  # save input file at LSsimulated
    lightning_dic, tdoa_warning = bckend.solver(lstrikes)           # solves tdoa geo/matrix
    if tdoa_warning:        
        print("Warning: tdoa_solve algorithm for more than 3 stations were not provided. Adopting geometrical solution (with 3 stations)")
    bckend.save(lightning_dic,filename=output_file_res)             # save results at LSresults
