import numpy as np
import datetime

# to access this module in command line:
# 
# $ python3
# >>> import sys
# >>> sys.path.append('/Users/robertmorgan/documents/git/rwm/circuits')
# >>> import circuit_analysis_tools as cat

# a circuit consists of nodes and elements.
# nodes have voltages
# elements have currents
# currents from start node to end node are positive
# 

############################
######## Node Class ########
############################

class Node:
    '''
    I want to do stuff with nodes too, so they should also exist
    ''' 
    def __init__(self, node_data, times):
        ''' 
        Initializes a node object.

        Inputs:

        node_data    dict, defines the nodes's ID and initial voltage, and
                     whether the node is grounded.
        '''

        # define the type and ID attributes
        self.type          = 'node'
        self.id            = node_data['id']
        self.is_ground     = node_data['is_ground']

        # create an array to log node voltage over time
        self.voltages_V    = np.full(len(times), np.nan)

        # set first voltage to the initial value
        self.voltages_V[0] = node_data['voltage_initial_V']

class Branch:
    '''
    a group of elements that are in series and thus have the same current
    '''
    def __init__(self, branch_data):
        '''
        initializes a branch object.
        '''
        # store branch parameters as attributes
        self.branch_name    = branch_data['branch_name']
        self.ref_element_id = branch_data['ref_element_id']

###############################
######## Element Class ########
###############################

class Element:
    def __init__(self, element_data, nodes, times):
        '''
        Initializes a circuit element object. 
        
        Inputs:

        element_data   dict, defnines the element's characteristics
        nodes          dict, the circuit's node objects. key is <node ID> 
                       and value is <node object>
        times          1d iterable, times at which circuit behavior
                       will be calculated. needed mostly for its length.
        '''

        #########################
        ### All element types ###
        #########################

        self.element_type  = element_data['type']
        self.element_id    = element_data['id']
        self.node_id_start = element_data['node_start']
        self.node_id_end   = element_data['node_end']    

        # assign start and end node objects as attributes
        self.node_start    = nodes[self.node_id_start]    
        self.node_end      = nodes[self.node_id_end]    

        # create an array to log element current over time
        self.currents_A      = np.full(len(times), np.nan)

        # set first current to the initial value
        self.currents_A[0]   = element_data['current_initial_A']

        ################
        ### Resistor ###
        ################
        if self.element_type == 'resistor':
            # store the resistance, in ohms
            self.resistance_ohm = element_data['resistance_ohm']


        ################
        ### Inductor ###
        ################
        elif self.element_type == 'inductor':
            # store the inductance, in henries
            self.inductance_H = element_data['inductance_H']


        #################
        ### Capacitor ###
        #################
        elif self.element_type == 'capacitor':
            # store the capacitance, in farads
            self.capacitance_F = element_data['capacitance_F']


        ###################################
        ### Source - Voltage - Constant ###
        ###################################
        elif self.element_type == 'source_voltage_constant':
            # store the source voltage, in volts
            self.voltage_V = element_data['voltage_V']


        #####################################
        ### Source - Voltage - Sinusoidal ###
        #####################################
        elif self.element_type == 'source_voltage_sine':
            # store the maximum (zero-to-peak) voltage
            self.vmax_V = element_data['voltage_max_V']
            # store the frequency of the sine source
            self.frequency_Hz = element_data['frequency_Hz']
            # calculate the voltage value as of the first time
            self.voltage_V = self.update_element(times[0])

        ############################
        ### Diode - Interpolated ###
        ############################
        elif self.element_type == 'diode_interpolated':
            # store the diode curve as an array
            self.diode_curve = np.array(element_data['diode_curve_[[V,A]]'])
            #print('Diode - Interpolated: self.diode_curve = ', 
            #      self.diode_curve)
            #print('Diode Voltages = ', self.diode_curve[:, 0])
            #print('Diode Currents = ', self.diode_curve[:, 1])
            self.diode_resistance = element_data['resistance_ohm']

        else:
            print('!!! WARNING: ELEMENT TYPE NOT RECOGNIZED !!!')

    def update_element(self, time_s):
        '''    
        Updates the voltage of a 'source_voltage_sine' object as a
        function of time
        '''
        if self.element_type in ['resistor', 'capacitor', 'inductor', 
                                 'source_voltage_const', 
                                 'diode_interpolated']:
            pass
        if self.element_type == 'diode_interpolated':
            pass
            
        elif self.element_type == 'source_voltage_sine':
            # Calculate and return the voltage as a function of time
            self.voltage_V = self.vmax_V * np.sin(2 * np.pi * time_s 
                                                  * self.frequency_Hz)


    def update_diode_interpolated(self):
        ''' 
        TO DO
        '''
        pass


###############################
######## Circuit Class ########
###############################

class Circuit:
    '''
    A circuit object represents the 
    network of circuit elements, and is used to coordinate interactions
    between the various circuit elements.
    '''
    def __init__(self, circuit_data, times):
        '''
        Initializes a circuit object. 
        
        The 'circuit_data' input contains all the information 
        we need in order to build the circuit network, and is loaded from a 
        json file defining the circuit. 

        The json file must include the
        following:
            circuit_elements    a list of dictionaries, in which each
                                dictionary describes a circuit element

        'times' is a 1d array of times at which to solve the circuit.

        '''
        
        ######################################
        ### store the inputs as attributes ###
        ######################################

        self.circuit_data = circuit_data
        self.times        = times
        self.times_log    = [] # for debugging purposes

        ### initialize an attribute to track the index and value of 'now'
        self.now_index = 0
        self.now_time  = times[0]


        #############################
        ### generate node objects ###
        #############################

        # create a list to store the node IDs and a dict to store the 
        # node objects
        self.node_ids = []
        self.nodes    = {}

        # generate a node object for each node defined in the JSON data
        for node_data in self.circuit_data['nodes']:
            # generate the node object and extract its node ID
            node    = Node(node_data, self.times)
            node_id = node.id

            # add the node ID to the list of node IDs
            self.node_ids.append(node_id)

            # add node object to the nodes dict, using ID as a key
            self.nodes[node_id] = node


        ################################
        ### generate element objects ### 
        ################################

        # create a list to store the element IDs and a dict to store the 
        # element objects
        self.element_ids = []
        self.elements    = {} 
        
        # generate an element object for each element defined in the JSON data
        for element_data in self.circuit_data['circuit_elements']:
            # generate the element object and extract its element ID
            element    = Element(element_data, self.nodes, self.times)
            element_id = element.element_id

            # add the element ID to the list of element IDs
            self.element_ids.append(element_id)

            # add the element object to the elements dict, using ID as a key
            self.elements[element_id] = element
        
            # assign the element its start node object and end node object
            # as element as attributes, for convenience
            element.node_start = self.nodes[element.node_id_start]
            element.node_end   = self.nodes[element.node_id_end]

        ###############################
        ### generate branch objects ###
        ###############################
        self.branch_names = []
        self.branches     = {}

        for branch_data in self.circuit_data['branches']:
            # generate the branch object and extract its branch ID
            branch      = Branch(branch_data)
            branch_name = branch.branch_name

             # add the element ID to the list of element IDs
            self.branch_names.append(branch.branch_name)

            # add the element object to the elements dict, using ID as a key
            self.branches[branch_name] = branch

        #######################################
        ### organize the elements and nodes ###
        #######################################

        # assemble the elements and nodes into a vector, so that we can
        # keep track of which node voltage and element current goes where
        # as we build and solve our matrix equation Ax=b, This vector x_key 
        # is constructed as a list of [element IDs, node IDs].
        # The order of element IDs that was passed to __init__ is maintained. 
        # The nodes are sorted alphabetically. 
        self.x_key = self.element_ids + self.node_ids


        ###################
        ### debug tools ###
        ###################

        
    def increment_time(self):
        '''
        a function to increment to the next time step by updating the 
        meaning of 'now', the index of the present time step, and the
        time elapsed since the previous time step.
        '''

        # increment index of current time by now
        self.now_index += 1 

        # update the meaning of 'now'
        self.now_time = self.times[self.now_index]

        # calculate the time step from previous time to now
        self.dt = self.times[self.now_index] - self.times[self.now_index-1]
        # WHAT HAPPENS IF T_0 != 0?


    def update_circuit(self):
        '''
        This function exists to update the behavior of nonlinear elements
        like diodes, and time-variant things like AC voltage sources, to
        establish the latest approximations and values for these elements
        as of the present time step. 
        '''
        # establish the latest value of "now" in case time has incremented
        t_now = self.now_time
        # update each element
        for key in self.element_ids:
            e = self.elements[key]
            e.update_element(t_now)


    def get_matrix(self):
        '''
        creates A and b in Ax=b based on the elements and nodes in self
        '''

        # reset A and b to zeros
        A_dim = len(self.x_key) # A has dimensions n by n; b has length n
        
        A = np.zeros((A_dim, A_dim))
        b = np.zeros(A_dim)

        #############################################################
        ### calculate and store the coefficients for each element ###
        #############################################################

        for i in range(len(self.elements)):
            e_id   = self.element_ids[i]
            e      = self.elements[e_id] # extract the circuit element object
            e_type = e.element_type
            
            # define coefficients c_I, c_Vsn, c_Ven, c_b such that
            #    c_I  * element_current 
            #  + c_ns * start_node_voltage
            #  + c_ne * end_node_voltage
            #  = c_b 

            ###################################
            ### Source - Voltage - Constant ###
            ###################################
            if e_type == 'source_voltage_constant':
                c_I  =  0
                c_ns = -1
                c_ne =  1
                c_b  =  e.voltage_V

            #####################################
            ### Source - Voltage - Sinusoidal ###
            #####################################
            elif e_type == 'source_voltage_sine':
                c_I  =  0
                c_ns = -1
                c_ne =  1
                c_b  =  e.voltage_V
 
            ################
            ### Resistor ###
            ################
            elif e.element_type == 'resistor':
                c_I  =  e.resistance_ohm
                c_ns =  -1
                c_ne =  1
                c_b  =  0

            ################
            ### Inductor ###
            ################
            elif e.element_type == 'inductor':
                # pull out some things that will be useful later

                # define a constant k_L = delta-t / (2 * L)
                k_L = self.dt / (2. * e.inductance_H)

                # calculate voltage from previous time step at start node
                # and end node
                t_index_prev      = self.now_index - 1
                V_start_node_prev = e.node_start.voltages_V[t_index_prev]
                V_end_node_prev   = e.node_end.voltages_V[t_index_prev]
                I_prev            = e.currents_A[t_index_prev] 
                
                # calculate the coefficicents
                c_I  =  1
                c_ns = -k_L
                c_ne =  k_L
                c_b  =  (
                    I_prev                     # previous current
                    + k_L * V_start_node_prev  # previous start voltage
                    - k_L * V_end_node_prev    # previous end voltage
                    )

            #################
            ### Capacitor ###
            #################

            elif e.element_type == 'capacitor':
                # pull out some things that will be useful later

                # define a constant k_C = delta-t / (2 * L)
                k_C = self.dt / (2. * e.capacitance_F)

                # calculate voltage from previous time step at start node
                # and end node
                t_index_prev      = self.now_index - 1
                V_start_node_prev = e.node_start.voltages_V[t_index_prev]
                V_end_node_prev   = e.node_end.voltages_V[t_index_prev] 
                I_prev            = e.currents_A[t_index_prev] 
                
                # calculate the coefficicents
                c_I  = -k_C
                c_ns =  1
                c_ne = -1
                c_b  =  (
                    k_C * I_prev         # previous current
                    + V_start_node_prev  # previous start voltage
                    - V_end_node_prev    # previous end voltage
                    )

            ############################
            ### Diode - Interpolated ###
            ############################

            elif e.element_type == 'diode_interpolated':
                # calculate voltage from previous time step at start node
                # and end node
                t_index_prev      = self.now_index - 1
                V_start_node_prev = e.node_start.voltages_V[t_index_prev]
                V_end_node_prev   = e.node_end.voltages_V[t_index_prev] 
                Vf_prev           = V_start_node_prev - V_end_node_prev

                # calculate the voltages from any prior iterations of
                # the present time step
                t_index_iter      = self.now_index
                V_start_node_iter = e.node_start.voltages_V[t_index_iter]
                V_end_node_iter   = e.node_end.voltages_V[t_index_iter] 
                Vf_iter           = V_start_node_iter - V_end_node_iter
                # print('V_start_node_iter = ', V_start_node_iter)
                # print('V_end_node_iter = ', V_end_node_iter)
                # print('Vf_iter = ', Vf_iter)

                
                # find the previous time step current based on the
                # previous time step Vf and the interpolated diode table
                I_prev_interp = np.interp(Vf_prev,
                                          e.diode_curve[:,0],
                                          e.diode_curve[:,1])

                # calculate the local diode curve slope. If available, use
                # any prior iterations of this time step. Otherrwise, use the
                # slope of the diode table at the prior time step voltage.

                # use prior iteration of this time step if available
                if not np.isnan(Vf_iter):
                    I_iter_interp = np.interp(Vf_iter,
                                              e.diode_curve[:,0],
                                              e.diode_curve[:,1])
                    delta_I = I_iter_interp - I_prev_interp
                    delta_V = Vf_iter - Vf_prev
                    if delta_V == 0:
                        k_d = 0
                    else:
                        k_d = delta_I / delta_V

                else:
                    # calculate the local slope using the previous time
                    # step voltage and the diode table
                    h    = 0.01 # (float) [V] increment for slope calculation
                    I_V  = np.interp(Vf_prev,
                                     e.diode_curve[:,0],
                                     e.diode_curve[:,1])

                    I_Vh = np.interp(Vf_prev + h,
                                     e.diode_curve[:,0],
                                     e.diode_curve[:,1])

                    k_d  = (I_Vh - I_V) / h
                    if k_d < 0:
                        print('k_d = ', k_d)

                R = e.diode_resistance

                # calculate the coefficients
                c_I  =  1
                c_ns = -k_d - (1. / R)
                c_ne =  k_d + (1. / R)
                c_b  =  I_prev_interp - (k_d * Vf_prev)


            ##############################
            ### Store the coefficients ###
            ##############################

            # update the elements in row i of A to reflect the coefficients
            # calculated above:
            A[i, self.x_key.index(e.element_id   )] = c_I
            A[i, self.x_key.index(e.node_id_start)] = c_ns
            A[i, self.x_key.index(e.node_id_end  )] = c_ne

            # update element i within the vector b
            b[i] = c_b



        ######################################################
        ### calculate and store the coefficients for nodes ###
        ######################################################

        # create a ground node tracker to ensure we only apply one
        ground_node_applied = False
        # note the last column index used for elements above
        j = len(self.elements) - 1

        # One node is not linearly independent, and should be omtited from
        # the system of equations. As a matter of convention, choose that 
        # node to be the last node in the list, except that if the last node
        # is the ground node, omit the second-to-last node from the system.

        #choose the set of node ids to include in the system
        node_ids_to_solve = self.node_ids[:]
        # remove the last node from the list, if it is not ground
        if self.nodes[node_ids_to_solve[-1]].is_ground == False:
            del node_ids_to_solve[-1]
        #otherwise, remove the second to last node
        else:
            del node_ids_to_solve[-2]

        # iterate through nodes, skipping the last node, and define their
        # A coefficients based on the KCL for that node. Note that b is
        # always zero, so no need to update that. Skip the last node, because
        # it is not linearly indpendent.
        # note that this means the ground node cannot be listed last

        for n_id in node_ids_to_solve:
            n    = self.nodes[n_id]   # extract the node object

            #     # increment to the next column
            j += 1

            
            ###############################################
            ### Construct the KCL equation for the node ###
            ###############################################
            # iterate through all the elements that could start or end here
            for e_id in self.elements:
                e = self.elements[e_id]
                # check whether this element starts at this node
                if e.node_id_start == n_id:
                    # if it does, set the 'arriving current' coeff to -1
                    A[j, self.x_key.index(e.element_id)] = -1
                # check whether this element ends at this node. note that
                # this is not "elif" to catch elements that start and end
                # at the same node. 
                if e.node_id_end == n_id: 
                    # if it does, set the 'arriving current' coeff to 1
                    A[j, self.x_key.index(e.element_id)] =  1

            ##############################################################
            ### Construct the Ground Node Equation if node is grounded ###
            ##############################################################
            if n.is_ground == True:
                if ground_node_applied == True:
                    # TO DO: raise an ereror, multiple grounds found
                    print('WARNING: MULTIPLE GROUNDS')  
                # make the correct row and column equal to 1, to set
                # up the equation 1 * V_ground_node = 0  
                A[len(self.x_key) - 1, self.x_key.index(n_id)] = 1
                ground_node_applied = True

            
        #print('A = ', A)
        #print('b = ', b)

        # update matrix attributes
        self.A = A
        self.b = b


    def solve(self):
        '''
        a function to solve for the current time step and store the results
        '''
        x = np.linalg.solve(self.A, self.b)
        
        # store currents and voltages to to element and node objects
        for i in range(len(x)):
            key = self.x_key[i]
            if key in self.element_ids:
                self.elements[key].currents_A[self.now_index] = x[i]
            elif key in self.node_ids:
                self.nodes[key].voltages_V[self.now_index] = x[i]

        # for debug
        self.x = x

############################################
######## Element Behavior Functions ########
############################################

# PORTED TO ELEMENT CLASS AS A METHOD
# def update_source_voltage_sine(element, time_s):
#     '''
#     TO DO: Where should this function live? should it be a method of any
#     class? 

#     TO DO: doing this here because it's a nice example of an update function
#     that can be called to update an element's behavior. But this one can
#     be pre-calculated based on times, so maybe it should be defined at
#     instantiation?

#     Updates the voltage of a 'source_voltage_sine' object as a
#     function of time

#     '''

#     # Calculate and return the voltage as a function of time
#     V_now = np.sin(2 * np.pi * time_s * element.frequency_Hz)
#     return V_now


def coefficients(element_or_node_object):
    '''
    input is an element or node object.
    '''
    pass

def get_timestamp():
    '''
    returns a string representing the date and time at which the function
    is called.

    output format is 'YYYY-MM-DD_HH-MM-SS'
    '''

    # create a raw date string with the format '2020-07-05 08:49:56.287390'
    datestring_raw = str(datetime.datetime.now())
    # now build a new version that omits the sub-second precision and changes
    # characters to be filename friendly
    datestring     = ''
    for i in range(datestring_raw.index('.')):
        if datestring_raw[i] == ' ':
            datestring += '_'
        elif datestring_raw[i] == ':':
            datestring += '-'
        else:
            datestring += datestring_raw[i]
    
    return datestring










