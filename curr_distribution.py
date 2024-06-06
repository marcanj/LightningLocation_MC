import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, LogFormatter
import json 


def curr_dist_parameters(data):
    """Fit a lognormal distribution to the data and return the parameters of the distribution"""
    data = np.abs(data)                                 # Convert 'data' to a numpy array if it's not already
    data = data.values.ravel()                          # Assuming we want a 1D array
    dataf = data[data >= 2e3]                           # filter data less than 2000 Amps
    shape, loc, scale = stats.lognorm.fit(dataf)        # Fit a lognormal distribution to your data
    return shape, loc, scale

def plot_histograms(data, param):
    """Plot the histogram of the data and the histogram of the synthetic data generated using the parameters"""
    plt.figure(figsize=(12, 6))
    data_generated = stats.lognorm.rvs(param['shape'], loc=param['loc'], scale=param['scale'], size=data.size)
    plt.hist(data_generated, bins=np.linspace(0, 200000, 100), alpha=0.5, density=True, label=f'{param["pol"]} synthetic data')
    plt.hist(data.values.ravel(), bins=np.linspace(0, 200000, 100), alpha=0.5, density=True, label='real data')
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.title('Comparison of Synthetic Data with Actual Data')
    plt.legend()
    plt.show()
    return data_generated

def plot_cumulative_distribution(data):
    """Plot the cumulative distribution of the data"""
    sorted_df = data.sort_values(by='peakcurrent', ascending=True)              # Sort the data
    cumulative_prob = np.arange(len(sorted_df), 0, -1) / len(sorted_df)         # Compute the cumulative distribution
    log_current = np.log10(sorted_df)                                           # Compute the log10 of the current
    
    # Plot
    plt.figure(figsize=(8, 6))
    plt.plot(sorted_df, cumulative_prob, marker='o', linestyle='-')
    plt.xlabel('log10(Current)')
    plt.xscale('log') 
    plt.ylabel('P(I > i)')
    plt.yscale('log') 
    plt.title('Cumulative Distribution of Current')
    plt.grid(True)
    plt.show()


"""
Main code: 
Fit a lognormal distribution to the data and plot the histograms of the data and the synthetic data generated using the parameters

"""
try:
    input_pos_strikes = 'posRS.csv'                         # single column csv file with peak current of positive strikes given in Amperes
    input_neg_strikes = 'negRS.csv'                         # single column csv file with peak current of negative strikes given in Amperes 
    output_file = 'curr_distribution_parameters.json'       # output file to store the parameters of the lognormal current distribution

    parameters = []

    data_pos = pd.read_csv(input_pos_strikes)
    shape, loc, scale = curr_dist_parameters(data_pos)
    parameters.append({'pol':'pos', 'shape':shape, 'loc':loc, 'scale':scale})

    data_neg = pd.read_csv(input_neg_strikes)
    shape, loc, scale = curr_dist_parameters(data_neg)
    parameters.append({'pol':'neg', 'shape':shape, 'loc':loc, 'scale':scale})

    with open(output_file, 'w') as json_file:
        json.dump(parameters, json_file)

    data_pos_gen = plot_histograms(data_pos, parameters[0])
    data_neg_gen = plot_histograms(-data_neg, parameters[1])

    #plot_cumulative_distribution(-data_neg)
    #plot_cumulative_distribution(data_pos)
except:
    print('Error: Check if the input files are correct and the data is in the correct format')
    print('Make sure the input files are in the same directory as this script, named posRS.csv and negRS.csv')
    print('Make sure the input files are single column csv files with peak current of positive and negative strikes given in Amperes')




