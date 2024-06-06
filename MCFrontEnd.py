from flask import Flask, render_template, request, redirect     # Import redirect
from MCBackend import MCBackend
from flask import Flask, session

# INPUT FILES 
input_stations_file = 'STA_network.json'                        # station network file
curr_dist_parameters_file = 'curr_distribution_parameters.json' # file with the parameters of the current distribution
svm_parameters_negLS = 'svm_peak_neg_300km.pkl'                # file with the parameters of the SVM model
svm_parameters_posLS = 'svm_peak_pos_300km.pkl'                # file with the parameters of the SVM model

# OUTPUT FILES
output_file_res = 'LSresults.json'                              # output file for results

app = Flask(__name__)
app.secret_key = b'VeryHardToGuessSecretKey123'                 # Secret key is not important to be revealed in this code

# initialize backend
bckend = MCBackend(input_stations_file, curr_dist_parameters_file, svm_parameters_negLS, svm_parameters_posLS)      

@app.route('/', methods=['GET', 'POST'])
def home():
    if request.method == 'POST':
        selected_region = request.form.get('region')            # Get the selected region from the form
        return redirect(f'/index?region={selected_region}')     # Redirect to the index page with the selected region
    regions = bckend.get_regions()
    return render_template('home.html', regions=regions)        # Render the home page with the list of regions

@app.route('/index')
def index():
    region = request.args.get('region', 'USA-FL')               # Get the selected region from the URL query parameters
    selected_stations = bckend.filter_stations_by_region(region) # Filter stations based on the selected region
    return render_template('index.html', available_stations=selected_stations) # Render the index page with the available stations

@app.route('/check', methods=['POST'])
def check():
    selected_stations = request.form.getlist('stations')
    selected_stations_info = [station for station in bckend.stations if station['name'] in selected_stations]
    avg_coordinates = bckend.calculate_average_coordinates(selected_stations_info)
    return render_template('map.html', stations=selected_stations_info, avg_coordinates=avg_coordinates)

@app.route('/calculate', methods=['POST'])
def calculate():
    data = request.json
    blue_stations = data.get('blue_stations', [])
    red_stations = data.get('red_stations', [])
    blue_stations_json = [{'lat': station['lat'], 'lon': station['lon']} for station in blue_stations]
    red_stations_json = [{'lat': station['lat'], 'lon': station['lon']} for station in red_stations]

    session['blue_sta'] = blue_stations_json                        # Store the blue and red stations data in the session
    session['red_sta'] = red_stations_json
    session['num_ls'] = int(data.get('num_ls', 100))                # default Num_LS = 100
    session['tick'] = data.get('tickSelected', True)                # for default, selects mode 'True' 

    return "Worked"

@app.route('/operation-finished')
def operation_finished():
    num_ls = session['num_ls']
    input_stations=session['red_sta']+session['blue_sta']
    bckend.mode = session['tick']                                   # True: operates ignoring SVM pattern. False: considers SVM detection pattern
    bckend.add_input_stations(input_stations)
    avg_coordinates = bckend.calculate_average_coordinates(input_stations)
    StrikeList = bckend.commander(num_ls)
    lightning_dic, tdoa_warning = bckend.solver(StrikeList)
    if tdoa_warning:
        warning_message = "Warning: tdoa_solve algorithm for more than 3 stations were not provided. Adopting geometrical solution (with 3 stations)"
    else:
        warning_message = None
    bckend.save(lightning_dic,filename='LSresults.json')
    return render_template('map_L.html', red_stations=session['red_sta'], blue_stations=session['blue_sta'], lightning=lightning_dic, avg_coordinates=avg_coordinates, warning_message=warning_message)

if __name__ == '__main__':
    app.run(debug=True)