# this is a sandbox to try building and simulating circuits

import numpy as np
import json
import circuit_analysis_tools as cat
import matplotlib.pyplot as plt
import matplotlib as mpl
import datetime

# Output preferences
save_results_to_file = False
save_figure_to_file  = False #currently requires save to file to also be true

# Time step parameters
tmax = 1 # end time [s]
dt = tmax/100000. # time step [s]
times = np.arange(-dt, tmax, dt)

# plot time interval
t_plot_start = 0.4
t_plot_end = 0.6

###################
### input files ###
###################

### json file with node and element data
#config_JSON_filename = 'example_circuit_11_diode.json'
config_JSON_filename = 'dynamo_circuit_1.json'

### circuit diagram image filaname, can be str or None
#circuit_diagram_filename = 'example_circuit_11_diode.png'
circuit_diagram_filename = 'dynamo_circuit_1.png'


################################################
### load the circuit data from the json file ###
################################################

with open(config_JSON_filename, 'r') as read_file:
    circuit_data = json.load(read_file)

#################################
### create the circuit object ###
#################################

mycircuit = cat.Circuit(circuit_data, times)

##################################
### iterate through time steps ###
##################################

for i in range(len(mycircuit.times) - 1):

    # we do not solve at the first time, instead we take it as the
    # initial conditions. so we first increment time by one step
    mycircuit.increment_time()

    # then we calculate new element behaviors based on time and 
    # each element's previous current and previous node voltages 
    mycircuit.update_circuit()

    # then we generate this timestep's matrix equations and solve them
    mycircuit.get_matrix()
    A = mycircuit.A
    A_time = mycircuit.now_time
    A_step = mycircuit.now_index
    mycircuit.solve()


    # option to use additional iterations without incrementing time
    # (sometimes helpful to avoid diode interpolation overshoot)
    additional_iterations = 0
    for i in range(additional_iterations):
        mycircuit.get_matrix()
        mycircuit.solve()

    
    # and then we do it all again for each remaining time step!    


####################################################################
### write the circuit data and simulation results to a JSON file ###
####################################################################

if save_results_to_file:
    # create dictionaries to hold results
    output_JSON_data = {}    # the staging dict that will be written to JSON file
    results = {}             # dict to hold all simulation results
    node_voltages_V = {}     # dict to hold voltage results
    element_currents_A = {}  # dict to hold element current results

    # Time
    # store the array of times to the results dict, converting the array
    # to a list so that it is JSON serializable
    results['times_s'] = mycircuit.times.tolist()

    # Nodes
    # compile node voltages and store to node results dict, converting
    # the voltages array to a list so that it is JSON serializable
    for n_id in mycircuit.node_ids:
        n_voltages = mycircuit.nodes[n_id].voltages_V.tolist()
        node_voltages_V[n_id] = n_voltages

    # store the dict of node voltages to the main results dict
    results['node_voltages_V'] = node_voltages_V

    # Elements
    # compile element currents and store to element results dict, converting
    # the currents array to a list so that it is JSON serializable
    for e_id in mycircuit.element_ids:
        e_currents = mycircuit.elements[e_id].currents_A.tolist()
        element_currents_A[e_id] = e_currents

    # store the dict of element currents to the main results dict
    results['element_currents_A'] = element_currents_A

    # put the circuit data and the results into the JSON staging dict
    output_JSON_data['circuit_config']     = mycircuit.circuit_data
    output_JSON_data['simulation_results'] = results

    # create a timestamp string to add to the file name
    filename_timestamp = cat.get_timestamp()
    output_filename = 'simulation_output/sim_output_' + filename_timestamp

    # write circuit data and and simulation results to a JSON file #

    with open(output_filename + '.json', 'w') as output_file:
        json.dump(output_JSON_data, output_file, indent=4)


#############################
### debugging and testing ###
#############################
if False: # PRINT STUFF FOR DEBUGGING
    #print('all element IDs:    ', (e.element_id for e in mycircuit.elements))
    print('element IDs: ', mycircuit.element_ids)
    print('node IDs:    ', mycircuit.node_ids)
    print('element and node IDs (x_key):  ', mycircuit.x_key)

    if 0:
        print('')
        print('planned times = ', mycircuit.times)
        print('logged times  = ', mycircuit.times_log)
        for key in mycircuit.x_key:
            print('values at %s = '%key, mycircuit.x_obj_from_id(key).data)

    print('\n')
    # print('L1 currents = ', mycircuit.x_obj_from_id('L1'))
    print('mycircuit.times = ', mycircuit.times)
    print('times =            ', times)


################
### Plotting ###
################

### Make a combined plot of all the elements and nodes ###
if True:
    # figsize format is (width_inches, height_inches)
    fig = plt.figure(figsize=(15,8), dpi=100, facecolor='white') # (18,10) for ext monitor
    ax_n  = fig.add_subplot(221)
    ax_e  = fig.add_subplot(223)
    ax_im = fig.add_subplot(122)
    #fig.set_size_inches(6.5, 5) 
    for n_id in mycircuit.node_ids:
        n_voltages = mycircuit.nodes[n_id].voltages_V
        ax_n.plot(mycircuit.times, n_voltages, label=n_id)


    # Plot element currents by branch
    for branch_name in mycircuit.branch_names:
        branch = mycircuit.branches[branch_name]
        e_id   = branch.ref_element_id
        e_currents = mycircuit.elements[e_id].currents_A
        ax_e.plot(mycircuit.times, e_currents, label=branch_name)

    # set up labels, titles, legends
    ax_n.set_xlabel('time [s]')
    ax_n.set_ylabel('voltage [V]')
    ax_n.set_title('Node Voltages')
    ax_n.legend(loc='upper right')
    ax_n.set_xlim([t_plot_start, t_plot_end])
    ax_e.set_xlabel('time [s]')
    ax_e.set_ylabel('current [A]')
    ax_e.set_title('Branch Currents')
    ax_e.legend(loc='upper right')
    ax_e.set_xlim([t_plot_start, t_plot_end])


    # Add circuit diagram image
    if circuit_diagram_filename != None:
        circuit_diagram = plt.imread(circuit_diagram_filename)
        ax_im.imshow(circuit_diagram)
        fig.text(.75, .99, 'Diagram Source: ' + circuit_diagram_filename, 
            fontsize='16', horizontalalignment='center', 
            verticalalignment='top')
    else:
        fig.text(.75, .5, 'no diagram provided', color='#707070',
            fontsize='32', horizontalalignment='center', 
            verticalalignment='center', rotation='45')
    ax_im.axis('off')
    fig.set_tight_layout(True)
    if save_figure_to_file:
        fig.savefig(output_filename + '.png', dpi=200)
    
### make a plot to inspect selected node voltages over time
if False:
    for n_id in ['f']:#['b', 'c', 'd', 'e', 'f', 'g']:
        n_voltages = mycircuit.nodes[n_id].voltages_V
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(mycircuit.times, n_voltages)
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Voltage [V]')
        ax.set_title('Voltage over Time for Node %s' %n_id)

### Make a separate current plot for each branch
if False:
    for branch_name in mycircuit.branch_names:
        branch = mycircuit.branches[branch_name]
        e_id   = branch.ref_element_id
        e_currents = mycircuit.elements[e_id].currents_A
        fig = plt.figure() 
        ax  = fig.add_subplot(111)
        ax.plot(mycircuit.times, e_currents)
        ax.set_xlabel('time [s]')
        ax.set_ylabel('current [A]')
        ax.set_title('Branch Current, [' + branch_name + ']')

### make a plot to inspect an element's V-I curve
if False:
    elements_to_plot = ['D1', 'D5', 'S1']
    fig = plt.figure(figsize=(6,9), dpi=100, facecolor='white')
    for i in range(len(elements_to_plot)):
        e_id = elements_to_plot[i]
        e_currents = mycircuit.elements[e_id].currents_A
        e_voltages = (mycircuit.elements[e_id].node_start.voltages_V 
                      - mycircuit.elements[e_id].node_end.voltages_V)
        #fig = plt.figure()
        ax_e = fig.add_subplot(len(elements_to_plot), 1, i + 1)
        ax_e.plot(e_voltages, e_currents)
        ax_e.set_xlabel('voltage [V]')
        ax_e.set_ylabel('current [A]')
        ax_e.set_title('Current vs Voltage for %s' %e_id)
        #ax_e.set_xlim([t_plot_start, t_plot_end])
    fig.set_tight_layout(True)

plt.show()
















