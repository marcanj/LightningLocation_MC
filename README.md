# LightningLocation_MC
The Monte Carlo method is a computational technique that involves using random sampling to simulate and analyse complex systems, allowing for the estimation of probabilities, uncertainties, and behaviour patterns through repeated iterations.

For assessing Lightning Location, the method involves simulating lightning strikes across a defined area, adding features to mimic real-world conditions. Detected strike locations from the station are compared to the simulated strikes, statistically evaluating system reliability.

The real lightning current distribution and station performance (detection pattern) are input modules of the tool and are modelled as described in the following documentation, allowing changes without interfering in the core code:

<aside>
⚠️ Input files:

</aside>

- **STA_network.json**: list of stations with name, region and coordinates;
- **svm_peak_pos_300km.pkl** and **svm_peak_pos_300km.pkl**: definition of values currentXdistance where detection happens;
- **curr_distribution_parameters.json**: lognormal fit parameters for generating current of random lightning strikes.
- **tdoa_solve.py**: algorithm to solve using the time difference of arrival technique, for four stations or more. (If not provided, the solution will consider only three stations and use the built in geometrical solution)

<aside>
⚠️ Output files:

</aside>

- **LSresults.json**: dictionary with all lightning strike data, timestamp, current, simulated location, stations that detected with respective ToA, stations coordinates, actual location, tdoa error.

# Py Codes

- **MCBackend**: All functions and mathematical logics used to solve lightning location, including geometrical and matrix solvers;
- **MCBackend_offline**: Commander for allowing execution of code without FrontEnd interface.
- **MCFrontEnd**: Code and folders that handle user-interface (jinja/javascript)


![Screenshot 2024-06-07 at 14 52 54](https://github.com/marcanj/LightningLocation_MC/assets/101256937/cbd128e3-ef51-4ceb-a205-72faa3563550)


## Current Distribution

Lightning strikes present a stochastic behaviour, as a result of their inherently random nature, influenced by atmospheric conditions and electrical properties. Similarly, the distribution of lightning currents varies stochastically due to factors like location, channel geometry, and atmospheric conductivity. For creating a model of lightning current distribution, more than 400,000 lightning strikes are separated by polarity, and their current distribution is modelled, using a lognormal fitting. That distribution is incorporated to the model as well as the probability of a single strike being positive or negative (approximately ~10% to be positive).

The code for this model receives as input two csv tables with a single column with lightning currents, a single file per polarity (“posRS.csv” and “negRS.csv”). The “curr_distribution.py” file process the data and can be altered independently of the Monte Carlo processing, as long as the output “curr_distribution_parameters.json” is provided.

### Positive and Negative Lightning Current Distribution
Have a look at negCurrDistribution.png and posCurrDistribution.png.


## Detection Pattern

Lightning generates electromagnetic fields that propagate over a range of distances, with the amplitude of these fields depending on the peak current of the lightning strike. The higher the peak current, the stronger the electromagnetic field. As electromagnetic waves propagate through space, they naturally experience attenuation, or decay, over distance, caused by factors such as absorption and scattering in the surrounding medium.
The detection station, on the other hand, has its own limiting factors for sensing, mainly antennas and amplifiers that play an important role in the final waveform delivered to software for the time of arrival determination. In order to contemplate that, a thunderstorm in a range within 300 km of the station with more than 5000 lightning strikes was used.
Using real data from a calibrating Lightning Locating System (LLS), and combining the performance of the detector and the peak-detector algorithm, the input data was classified as ‘detected’ or ‘not detected’, and through a scatter plot current X distance, it was possible to observe the formation of detected/non-detected regions, based on the currents and distances of the signals.

Have a look at svm_peak_pos/neg_300km.png files.
