# -*- coding: utf-8 -*-

#%% import packages

from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.ticker as ticker

from scipy.integrate import quad

import pandas as pd

#%% Chemical class
class Chemical:
    def __init__(self, name, concentration=0):
        self.name = name
        self.concentration = concentration

    def set_concentration(self, concentration):
        self.concentration = concentration

    def get_concentration(self):
        return self.concentration

#%% Tank class
class Tank:
    def __init__(self, max_volume, current_volume=0, chemicals=None, source=None, dest=None):
        self.max_volume = max_volume
        self.current_volume = current_volume
        self.chemicals = chemicals if chemicals is not None else []
        # limited to one source and one destination
        self.source = source if chemicals is not None else []
        self.dest = dest if chemicals is not None else []
        
    def get_volume(self):
        return self.current_volume
    
    def set_volume(self, new_volume):
        self.current_volume = new_volume
    
    def get_chemicals(self):
        return self.chemicals
    
    def add_chemicals(self, chemical):
        return self.chemicals.append(chemical)
    
    def remove_chemicals(self, chemical):
        return self.chemicals.remove(chemical)
    
    def set_concentration(self, chemical_name, new_concentration):
        for chemical in self.chemicals:
            if chemical.name == chemical_name:
                chemical.set_concentration(new_concentration)
                break
    
    def get_concentrations(self):
        concentrations = {chem.name: chem.get_concentration() for chem in self.chemicals}
        return concentrations
    
    def get_masses(self):
        masses = {chem.name: chem.get_concentration() * self.current_volume for chem in self.chemicals}
        return masses
    
    def set_source(self, source):
        self.source = source

    def get_source(self):
        return self.source
    
    def set_dest(self, dest):
        self.dest = dest

    def get_dest(self):
        return self.dest
        
#%% Adsorber class
class Adsorber:
    def __init__(self, name, max_volume, current_volume=0, chemicals=None, C_bulk=None, C_bound=None, source=None, dest=None, doneRRtime=None):
        self.name = name
        self.max_volume = max_volume
        self.current_volume = current_volume
        self.chemicals = chemicals if chemicals is not None else []
        self.C_bulk = C_bulk if C_bulk is not None else np.zeros(K)  # Assuming K compartments
        self.C_bound = C_bound if C_bound is not None else np.zeros(K)  # Assuming K compartments
        self.source = source
        self.dest = dest
        self.doneRRtime = doneRRtime if doneRRtime is not None else 0 # Set to 0 unles given
        self.inRR = False
        
    def get_volume(self):
        return self.current_volume
    
    def set_volume(self, new_volume):
        self.current_volume = new_volume
    
    def get_chemicals(self):
        return self.chemicals
    
    def add_chemicals(self, chemical):
        return self.chemicals.append(chemical)
    
    def remove_chemicals(self, chemical):
        return self.chemicals.remove(chemical)
    
    def set_concentration(self, chemical_name, new_concentration):
        for chemical in self.chemicals:
            if chemical.name == chemical_name:
                chemical.set_concentration(new_concentration)
                break
    
    def get_concentrations(self):
        concentrations = {chem.name: chem.get_concentration() for chem in self.chemicals}
        return concentrations
    
    def get_masses(self):
        masses = {chem.name: chem.get_concentration() * self.current_volume for chem in self.chemicals}
        return masses
    
    def set_source(self, source):
        self.source = source

    def get_source(self):
        return self.source
    
    def set_dest(self, dest):
        self.dest = dest

    def get_dest(self):
        return self.dest
    
    def is_full(self, feed_conc, threshold):
        """ 
        Check if the adsorber is full based on the bulk concentration 
        of the last cell being greater than a specified threshold of the feedstock conc.
        """
        # Assuming C_bound is an array representing bound concentration in each cell
        # and K is the total number of cells
        return self.C_bulk[-1] >= threshold * feed_conc  # Check the last cell
    
    def set_doneRRtime(self, doneRRtime):
        self.doneRRtime = doneRRtime

    def get_RRtime(self):
        return self.doneRRtime
    
    def reset_states(self, current_time):
        # print(f"Resetting states for {self.name} at time {current_time} min.")
        self.C_bulk.fill(0)  
        self.C_bound.fill(0)  
        for chemical in self.chemicals:
            chemical.set_concentration(0)  
        # print(f"States after reset for {self.name}: C_bulk={self.C_bulk}, C_bound={self.C_bound}")
        self.doneRRtime = 0
        self.inRR = False

#%% Pump class
class Pump:
    def __init__(self, source=None, dest=None, flow_rate=0):
        # limited to one source and one destination
        self.source = source
        self.dest = dest
        self.flow_rate = flow_rate
        
    def set_source(self, source):
        self.source = source

    def get_source(self):
        return self.source

    def set_flow_rate(self, flow_rate):
        self.flow_rate = flow_rate

    def get_flow_rate(self):
        if isinstance(self.source,Valve) and self.source.is_open == False:
            return 0
        if isinstance(self.dest,Valve) and not self.dest.is_open:
            return 0
        if isinstance(self.source,Tank) and self.source.current_volume < 1e-6:
            return 0
        return self.flow_rate

    def set_dest(self, dest):
        self.dest = dest

    def get_dest(self):
        return self.dest
    
    def get_concentrations(self):
        return self.source.get_concentrations()
    
    def get_chemicals(self):
        return self.source.get_chemicals()

#%% Valve class
class Valve:
    def __init__(self, source=None, dest=None):
        self.is_open = True
        # limited to one source and one destination
        self.source = source
        self.dest = dest

    def open_valve(self):
        self.is_open = True

    def close_valve(self):
        self.is_open = False

    def is_valve_open(self):
        return self.is_open

    def connect(self, source, dest):
        self.source = source
        self.dest = dest
        
    def set_source(self, source):
        self.source = source

    def get_source(self):
        return self.source

    def set_dest(self, dest):
        self.dest = dest

    def get_dest(self):
        return self.dest

    def get_flow_rate(self):
        if self.is_open == True:
            return self.source.get_flow_rate()
        else: 
            return 0
    
    def get_concentrations(self):
        return self.source.get_concentrations()
    
    def get_chemicals(self):
        return self.source.get_chemicals()
    
#%% MixingPoint class
class MixingPoint:
    def __init__(self):
        self.sources = []  # List to store sources (e.g., valves)
        self.destinations = []  # List to store destinations (e.g., valves)
    
    def add_source(self, source):
        self.sources.append(source)
    
    def remove_source(self, source):
        if source in self.sources:
            self.sources.remove(source)
    
    def add_destination(self, destination):
        self.destinations.append(destination)
    
    def remove_destination(self, destination):
        if destination in self.destinations:
            self.destinations.remove(destination)
    
    def get_flow_rate(self):
        # Count destinations with non-zero flow rate
        active_destinations = [dest for dest in self.destinations if dest.get_flow_rate() > 0]
        num_active_destinations = len(active_destinations)
        
        if num_active_destinations == 0:
            print("No active destinations")
            return None
        
        total_flow_rate_in = sum(source.get_flow_rate() for source in self.sources)
        
        # Assuming fluid resistance is equal to all active outlets
        return total_flow_rate_in / num_active_destinations if total_flow_rate_in > 0 else 0
    
    def get_concentrations(self):
        active_destinations = [dest for dest in self.destinations if dest.get_flow_rate() > 0]
        num_active_destinations = len(active_destinations)
        total_flow_rate_in = num_active_destinations * self.get_flow_rate()
        combined_data = {}
        
        for source in self.sources:
            for chemical, concentration in source.get_concentrations().items():
                if chemical not in combined_data:
                    combined_data[chemical] = 0.0
                if total_flow_rate_in != 0:
                    combined_data[chemical] += (concentration * source.get_flow_rate() / total_flow_rate_in)
                else:
                    combined_data[chemical] = 0.0
        
        return combined_data
    
#%% tank_system_ode function

def tank_system_ode(y, t, tanks):
    
    # if t > 5:
    #     valve2.close_valve()
        
    # if t > 10:
    #     valve1.close_valve()
    
    # Create a list to store the derivatives for all state variables (concentrations and volumes)
    dydt = []

    num_chemicals = max(len(tank.get_chemicals()) for tank in tanks)

    for i, tank in enumerate(tanks):
        tank_vol = y[i * (1 + num_chemicals) + num_chemicals]
        # Calculate the change in volume for the tank
        vol_in_i = 0
        vol_out_i = 0

        if tank.source is not None:
            vol_in_i += tank.source.get_flow_rate()

        if tank.dest is not None:
            vol_out_i += tank.dest.get_flow_rate()

        # Calculate the volume change for the tank
        dVdt = vol_in_i - vol_out_i
        
        # Enforce volume limits (current_volume <= max_volume)
        if tank_vol > tank.max_volume:
            dVdt = 0

        # Calculate the change in concentration for each chemical in the tank
        for j, chemical in enumerate(tank.get_chemicals()):
            tank_chemical_conc = y[i * (1 + num_chemicals) + j]
            
            mass_in_j = 0  # mass flow in to tank_i
            mass_out_j = 0  # mass flow out of tank_i
                        
            if tank.source is not None:
                # Calculate the total source flow rates for this chemical
                if tank.source.get_concentrations().get(chemical.name):
                    mass_in_j += tank.source.get_flow_rate() * tank.source.get_concentrations().get(chemical.name)
                    # print('Source C: ', tank.source.get_concentrations().get(chemical.name))

            if tank.dest is not None:
                # Calculate the total destination flow rates for this chemical
                if tank.get_concentrations().get(chemical.name):
                    mass_out_j += tank.dest.get_flow_rate() * tank_chemical_conc
                    # print('Tank/exit C: ', tank.get_concentrations().get(chemical.name))

            # Calculate the concentration change for the chemical in this tank
            if tank_vol != 0:
                dCdt_chemical = (1 / tank_vol) * (mass_in_j - mass_out_j - tank_chemical_conc * dVdt)
                if tank_vol > tank.max_volume:
                    # If full, mass in - mass out - mass overflowing --> mass_overflowing = (flow_in - flow_out)*C_tank
                    dCdt_chemical = (1 / tank_vol) * (mass_in_j - mass_out_j - (tank.source.get_flow_rate() - tank.dest.get_flow_rate())*tank_chemical_conc)
            else:
                dCdt_chemical = 0  # Set a default value when concentration is zero
            
            dydt.append(dCdt_chemical)

        dydt.append(dVdt)

    return dydt

#%% adsorber capacity function

def C_capacity_func(max_cycles,switch_times,adsorbers):
    """
    Calculate the column capacity as a function of uses.
    You can modify this function based on the specific behavior of your system.
    """

    num_adsorbers = len(adsorbers)
    adsorber_uses = int(switch_times/num_adsorbers)
    return 209.82*(1-adsorber_uses/max_cycles)


#%% establish and solve ADSORBER simulation

K = 10
porosity = 0.32

k_a = 0.104
k_d = 0.022


#%% initial condition function

def create_initial_conditions_vector(adsorbers, K):
    """
    Creates a vector for the initial conditions of the adsorber system ODE.

    Args:
    adsorbers (list): List of Adsorber objects.
    K (int): Number of compartments in each adsorber.

    Returns:
    numpy.ndarray: Array of initial conditions for the ODE system.
    """
    initial_conditions = []
    
    for adsorber in adsorbers:
        for k in range(K):
            initial_conditions.append(adsorber.C_bulk[k])  # C_bulk for compartment k
            initial_conditions.append(adsorber.C_bound[k]) # C_bound for compartment k

    return np.array(initial_conditions)

# Tank initial conditions (tank1, which is the only one integrated)
def tank_initial_conditions_vector(tanks, tank_chem):
    initial_conditions = []

    for tank in tanks:
        initial_conditions.append(tank_chem.get_concentration()) # Concentration of tank 1
        initial_conditions.append(tank.get_volume()) # Volume of tank 1

    return np.array(initial_conditions)

# Example usage:
# conc_ads_bulk_bound_0_test = create_initial_conditions_vector(adsorbers, K)

#%% unzip odeint results per adsorber and bulk/bound

def update_ads_bulk_bound(conc_ads_bulk_bound, adsorbers, K):
    """
    Unzips the conc_ads_bulk_bound vector to update the state of the adsorbers per the adsorber list.

    Args:
    conc_ads_bulk_bound (numpy.ndarray): The flattened array containing C_bulk and C_bound values.
    adsorbers (list): List of adsorbers.
    K (int): Number of compartments in each adsorber.
    """
    # num_adsorbers = len(adsorbers)
    latest_sim_result = conc_ads_bulk_bound[-1,:]
    
    for ads_ind,adsorber in enumerate(adsorbers):
        start_idx = ads_ind * K * 2
        end_idx = start_idx + K * 2
        
        adsorber_bulk_bound = latest_sim_result[start_idx:end_idx]
        
        adsorber_bulk = np.zeros(K)
        adsorber_bound = np.zeros(K)
        
        for k in range(K):
            adsorber_bulk[k] = adsorber_bulk_bound[k*2]
            adsorber_bound[k] = adsorber_bulk_bound[1+k*2]
            
        adsorber.C_bulk = adsorber_bulk
        adsorber.C_bound = adsorber_bound
        
        adsorber.chemicals[0].set_concentration(adsorber_bulk[-1])
        
# update_ads_bulk_bound(conc_ads_bulk_bound, adsorbers, K)


#%% check if there is waste and calculate the amount for the time interval
def waste_calc(time1, time2, conc1, conc2, waste_list):
    # determine if the column receiving loaded column's breakthrough also has breakthrough
    # use trapezoid rule to approximate area under curve
    waste_amount = (time2-time1)*(conc1+conc2)*0.5
    waste_list.append(waste_amount)

def mass_calc(time1, time2, conc1, conc2, mass_list):
    # determine the mass fed over the time interval given
    # use trapezoid rule for area under inlet concentration curve
    mass_in = (time2-time1)*(conc1+conc2)*0.5
    mass_list.append(mass_in)


# set up the functions that require concentration (which will be varied)
def inlet_feed_concentration(t, conc, A):
    """
    Calculate the inlet feed concentration as a function of time.
    You can modify this function based on the specific behavior of your system.
    """
    # #Example: a simple linear increase followed by a constant value
    return A * math.sin(math.pi / 200 * (t)) + conc
    

def inlet_flow_function(t, nom_flow, A):
    """
    Calculate the inlet flow rate as a function of time.
    Modify based on specific system behavior.
    """
    return A * math.sin(math.pi / 200 * (t)) + nom_flow

def adsorber_system_ode(y, t, adsorbers, tank1C):
                        
    chemical_A_tank1.concentration = tank1C
    
    adsorber_volume = 1 # mL
    
    V_bulk = adsorber_volume*porosity/K # mL
    V_resin = adsorber_volume*(1-porosity)/K    
    
    # Create a list to store the state variable derivatives for each adsorber in list adsorbers
    # (states are bulk & bound concentrations in K cells)
    dydt = [] # list for dydt for each adsorber

    # num_chemicals = max(len(tank.get_chemicals()) for tank in tanks)

    for i, adsorber in enumerate(adsorbers):
        
        # two groups: C_A^bulk,k & C_A^bound,k
        y_groups = np.reshape(y[i*K*2:(i*K*2+K*2)], (K,2))
        dydt_groups = np.zeros_like(y_groups)
        
        for k in range(K):
            
            q = 0
            if adsorber.source is not None:
                q += adsorber.source.get_flow_rate()
                bulk_conc_in = adsorber.source.get_concentrations().get('A')
                
            bulk_conc_k = y_groups[k,0]
            bound_conc_k = y_groups[k,1]
            bulk_conc_k_minus = y_groups[k-1,0]
            
            if k == 0:
                
                dydt_groups[k,0] = ((q*(bulk_conc_in - bulk_conc_k)) / V_bulk
                                    - (k_a * bulk_conc_k * (C_capacity - bound_conc_k) - k_d * bound_conc_k) * V_resin/V_bulk)
            
            else:
            
                dydt_groups[k,0] = ((q*(bulk_conc_k_minus - bulk_conc_k)) / V_bulk
                                    - (k_a * bulk_conc_k * (C_capacity - bound_conc_k) - k_d * bound_conc_k) * V_resin/V_bulk)
                
            # EQUATIONS FOR C^BOUND
            
            dydt_groups[k,1] = (k_a * bulk_conc_k * (C_capacity - bound_conc_k) - k_d * bound_conc_k)
            
            # Update the adsorber's C_bulk and C_bound
            adsorber.C_bulk[k] = y_groups[k, 0]
            adsorber.C_bound[k] = y_groups[k, 1]
        
        dydt.append(dydt_groups.flatten())

    return np.array(dydt).flatten()

# Initialize separate chemical objects for each tank

chemical_A_tank0 = Chemical(name="A", concentration=2.5)
chemical_A_tank1 = Chemical(name="A", concentration=2.5)
chemical_A_adsorber1 = Chemical(name="A",concentration=0)
chemical_A_adsorber2 = Chemical(name="A",concentration=0)
chemical_A_adsorber3 = Chemical(name="A",concentration=0)
chemical_A_adsorber4 = Chemical(name="A",concentration=0)

# Create the Tanks with Chemicals
# NOTE - current_volume CANNOT BE ZERO!
tank0 = Tank(max_volume=5000, current_volume=5000, chemicals=[chemical_A_tank0], source=None, dest=None)
tank1 = Tank(max_volume=250, current_volume=150, chemicals=[chemical_A_tank1], source=None, dest=None)

# Create the Valves
valve0 = Valve()
valve1 = Valve()
valve2 = Valve()
valve3 = Valve()
valve4 = Valve()

tanks = [tank1]
tanks0 = [tank0]
valves = [valve1,valve2]

# ------------------------------------------------ #


#%% Run the simulation with the system above



# Set up the conditions for the scenario to run
flows = [0.5, 1, 1.5, 2]
concs = [2.5, 2.75, 3]
threshs = [0.7]
flow_amp = 0.2 
conc_amp = -0.2    
RRtime = 40
num_cols = 4

# Binding capacity
Ks = 100 # mg/mL
col_vols = np.zeros((len(flows), len(concs)))

# Initialize arrays for stored variables
wastes = np.zeros((len(flows), len(concs), len(threshs)))
overflow_wastes = np.zeros((len(flows), len(concs), len(threshs)))
total_waste = np.zeros((len(flows), len(concs), len(threshs)))
percent_overflow_waste = np.zeros((len(flows), len(concs), len(threshs)))
waste_percent = np.zeros((len(flows), len(concs), len(threshs)))
prods = np.zeros((len(flows), len(concs), len(threshs)))
masses = np.zeros((len(flows), len(concs), len(threshs)))
total_masses = np.zeros((len(flows), len(concs), len(threshs)))
runs = 0

# Cycle through the flow rates input
for flowind, flow in enumerate(flows):
    
    # Create the Pumps, which rely on the flowrate input
    # Currently, always tank/adsorber -> pump -> valve
    pump0 = Pump(flow_rate=flow)
    pump1 = Pump(flow_rate=flow)
    pump2 = Pump(flow_rate=flow)
    pump3 = Pump(flow_rate=flow)
    pump4 = Pump(flow_rate=flow)
    pumps = [pump1,pump2]

    # Set up configuration for the system
    tank0.dest = pump0
    pump0.source = tank0
    pump0.dest = valve0
    valve0.source = pump0

    valve0.dest = tank1
    tank1.source = valve0

    tank1.dest = pump1
    pump1.source = tank1
    pump1.dest = valve1
    valve1.source = pump1

    pump2.dest = valve2
    valve2.source = pump2

    pump3.dest = valve3
    valve3.source = pump3

    pump4.dest = valve4
    valve4.source = pump4

    # the sequence is pump -> valve -> adsorber
    # the sequence starts with tank1, then goes to the first pump-valve-adsorber
    # after the first pump-valve-adsorber sequence, it goes to the next per the adsorbers list



    # Cycle through the concentration value inputs
    for concind, conc in enumerate(concs):
            # Use flow rate to determine column size
            col_vols[flowind, concind] = RRtime*conc*flow/Ks
            print('Column volume: %.2f' % col_vols[flowind, concind])
            
            KstoC = Ks/conc
            print('Ks to Cf ratio: %.2f' % KstoC)

            # Create the Adsorbers with Chemicals
            adsorber1 = Adsorber(name="Ads1",max_volume=col_vols[flowind, concind],current_volume=col_vols[flowind, concind],chemicals=[chemical_A_adsorber1],source=None,dest=None)
            adsorber2 = Adsorber(name="Ads2",max_volume=col_vols[flowind, concind],current_volume=col_vols[flowind, concind],chemicals=[chemical_A_adsorber2],source=None,dest=None)
            adsorber3 = Adsorber(name="Ads3",max_volume=col_vols[flowind, concind],current_volume=col_vols[flowind, concind],chemicals=[chemical_A_adsorber3],source=None,dest=None)
            adsorber4 = Adsorber(name="Ads4",max_volume=col_vols[flowind, concind],current_volume=col_vols[flowind, concind],chemicals=[chemical_A_adsorber4],source=None,dest=None)
            adsorbers = [adsorber1,adsorber2,adsorber3,adsorber4]

            valve1.dest = adsorbers[0]
            adsorbers[0].source = valve1
            adsorbers[0].dest = pump2
            pump2.source = adsorbers[0]

            valve2.dest = adsorbers[1]
            adsorbers[1].source = valve2
            adsorbers[1].dest = pump3
            pump3.source = adsorbers[1]

            valve3.dest = adsorbers[2]
            adsorbers[2].source = valve3
            adsorbers[2].dest = pump4
            pump4.source = adsorbers[2]

            valve4.dest = adsorbers[3]
            adsorbers[3].source = valve4


            for thresind, thresh in enumerate(threshs):

                # Reset the columns to avoid starting with partially filled adsorbers, cycling through each
                for adsorber in adsorbers:
                    adsorber.reset_states(0) # reset each adsorber at a time of 0

                # Reassign capacity to starting value for each trial
                C_capacity = 209.82             
                y0 = np.zeros((K, 2*len(adsorbers))).flatten()

                # Set the initial tank concentrations based on input values
                tank0.set_concentration("A", conc)
                tank1.set_concentration("A", conc)

                #%% Run the simulation with the system above
                time_steps = 60*8
                switch_thresh = thresh
                switch_duration = 0

                # Define time points for integration
                t_array = np.linspace(0, time_steps, (time_steps+1))
                delta_t = t_array[1] - t_array[0]
                    
                # Create initial concentrations for odeint
                conc_vol_tank0_0 = tank_initial_conditions_vector(tanks0, chemical_A_tank0)
                conc_vol_tank1_0 = tank_initial_conditions_vector(tanks, chemical_A_tank1)
                conc_ads_bulk_bound_0 = create_initial_conditions_vector(adsorbers, K)

                sim_results = np.zeros(shape=(len(t_array),len(adsorbers)+1))

                adsorber_to_column = {
                    'Ads1': 1, # 'Ads1' was the assigned adsorber name when the adsorber was defined above
                    'Ads2': 2, # Note that column 0 represents the time, so start with 'Ads1' : 1, not 0
                    'Ads3': 3,
                    'Ads4': 4
                }

                # Make vectors for the adsorbers to track loading, breakthrough, regeneration
                # Regeneration = -1
                # Waiting = 0
                # Breakthrough receiving = 1
                # Loading = 2
                Ads_stage_array = np.zeros(shape=(len(t_array), len(adsorbers)))
                
                # Initialize empty lists for storing values
                switch_times = []
                col_switches = 0

                waste_amounts = []
                overflow_waste_amounts = []

                mass_amounts = []
                total_mass_amounts = []
                
                tank1concs = []
                tank1vols = []
                tank0vols = []
                tank0concs = []

                # Add the starting volume and concentration of Tank1 and Tank0 to the lists
                tank1vols.append(conc_vol_tank1_0[1])
                tank1concs.append(conc_vol_tank1_0[0])

                tank0vols.append(conc_vol_tank0_0[1])
                tank0concs.append(conc_vol_tank0_0[0])

                # Perform time-stepping loop
                for t_ind, t in enumerate(t_array):
                    # Set the flow rates to flow(t)
                    flow_t = inlet_flow_function(t, flow, flow_amp)
                    if flow_t <= 0:
                        print(f"WANRING: FLOW RATE FROM BIOREACTOR <= 0!\n at t = {t}, Q = {flow}, C = {conc}, T = {thresh}")
                    pump0.flow_rate = flow_t
                    pump1.flow_rate = flow_t
                    pump2.flow_rate = flow_t
                    pump3.flow_rate = flow_t

                    # Set Tank0 concentration to match the inlet concentration changing with time
                    tank0.set_concentration("A", inlet_feed_concentration(t, conc, conc_amp))
                    conc_vol_tank0_0 = [inlet_feed_concentration(t, conc, conc_amp), tank0vols[-1]]
                                    
                    # Integrate the ODE system for Tank0 for this time step
                    odeint_time = [t, (t + delta_t)]
                    # Result is 2 x 2 with first row initial conditions and second row conditions at t + dt as [concentration, volume]
                    tank0_soln = odeint(tank_system_ode, conc_vol_tank0_0, odeint_time, args=(tanks0, ))
                    # Store the results for Tank0 integration
                    tank0concs.append(tank0_soln[-1,0])
                    tank0vols.append(tank0_soln[-1,1])

                    # Determine total mass entering surge vessel/leaving tank0
                    mass_calc(t, t+delta_t, tank0concs[-1], tank0concs[-2], total_mass_amounts)

                    
                    # Check what adsorbers are in RR, and if they are, check if they have completed at this t
                    for check_ind, ads_check in enumerate(adsorbers):

                        # Get index in the adsorber stage array for the adsorber
                        ads_stage_ind = adsorber_to_column[ads_check.name] - 1

                        # See if adsorber is in RR
                        if ads_check.inRR == True:
                            # Check if it has completed RR
                            if t >= ads_check.doneRRtime:
                                ads_check.inRR = False
                                # Set the time to be done RR to 0
                                ads_check.set_doneRRtime(0)
                                # Leave as 0 in adsorber stage array and open the valve
                                ads_check.source.open_valve()
                            else:
                                # Adsorber not done => add RR code to array, make source flow 0 (close valve)
                                ads_check.source.close_valve()
                                Ads_stage_array[t_ind, ads_stage_ind] = -1
                    

                    # Now, check if first adsorber is/still is in RR or if it can be loaded
                    if adsorbers[0].inRR == True:
                        # Stop any flow out of tank1 if not ready to load
                        tank1.dest.flow_rate = 0
                                                
                    else:
                        # Set flow to 'flow' if it can be loaded
                        tank1.dest.flow_rate = flow_t
                        
                    

                    # Check that there is enough volume left in the surge vessel to run
                    if tank1vols[-1] <= flow_t*1:
                        # Stop tank from drying up: set the flow out of the tank to 0
                        tank1.dest.flow_rate = 0
                        # print('Surge vessel outlet closed to prevent emptying at %i mins' % t)
                    else:
                        tank1.dest.flow_rate = flow_t


                    # Use the resulting flow rates to solve Tank1 ODE
                    tank1_soln = odeint(tank_system_ode, conc_vol_tank1_0, odeint_time, args=(tanks, ))
                    # Add solved values to lists of the concentrations and volumes of Tank1 and Tank0
                    tank1concs.append(tank1_soln[-1,0])
                    tank1vols.append(tank1_soln[-1,1])

                    # Solve for any waste from overflow
                    if tank1vols[-1] >= tank1.max_volume:
                        intgl_spill = []
                        # Integrate from t to t+delta t using tank 1 concentrations
                        waste_calc(t, t+delta_t, tank1concs[-1], tank1concs[-2], intgl_spill)
                        # Multiply by the flow rate of spill
                        overflow_waste_amounts.append((tank1.source.get_flow_rate() - tank1.dest.get_flow_rate())*intgl_spill[-1])
                
                    # Check if second adsorber can receive breakthrough
                    if adsorbers[1].inRR == True:
                        # Stop any flow out of first adsorber if not ready to load
                        adsorbers[1].source.close_valve()
                    else:
                        # Set flow rate to flow 
                        adsorbers[1].source.open_valve()

                    # Close valves for adsorbers[2] and adsorbers[3] regardless
                    adsorbers[2].source.close_valve()
                    adsorbers[3].source.close_valve()

                    # Solve the appropriate adsorbers based on which are being loaded
                    conc_ads_bulk_bound = odeint(adsorber_system_ode, conc_ads_bulk_bound_0, odeint_time, args=(adsorbers, tank1_soln[-1,0],))
                    update_ads_bulk_bound(conc_ads_bulk_bound, adsorbers, K)
                                                                             
                    # Record results in sim_results
                    sim_results[t_ind,0] = t
                                            
                    # Record results in sim_results using the dictionary
                    for adsorber in adsorbers:
                        column_index = adsorber_to_column[adsorber.name]
                        sim_results[t_ind, column_index] = adsorber.C_bulk[-1]
                    


                    # If there is any loading
                    if adsorbers[0].inRR == False:
                        load_arr_ind = adsorber_to_column[adsorbers[0].name] - 1
                        # Set the stage array value to 2 for loading
                        Ads_stage_array[t_ind, load_arr_ind] = 2
                        # Integrate for mass fed only if adsorbers[0] is being loaded (using tank1 concentration reaching adsorbers)
                        mass_calc(t, t+delta_t, tank1concs[-2], tank1concs[-1], mass_amounts)
                            
                        # SWITCH TEST - Check if adsorber is full
                        if adsorbers[0].is_full(feed_conc = conc, threshold=switch_thresh):
                            # Use dictionary to retrieve column index in sim_results for the breakthrough receiving adsorber
                            if adsorbers[1].inRR == False:
                                # If it is available, breakthrough is received by adsorbers[1]
                                bkt_col_index = adsorber_to_column[adsorbers[1].name]
                                # Set the stage array to 1 for receiving breakthrough
                                Ads_stage_array[t_ind, bkt_col_index - 1] = 1
                                                        
                            elif adsorbers[1].inRR == True:
                                # If not available, the product exiting adsorbers[0] is waste
                                bkt_col_index = adsorber_to_column[adsorbers[0].name]

                            # Do integral for waste on the adsorber
                            if sim_results[t_ind, bkt_col_index] >= 0.00001:
                                waste_calc(t, t+delta_t, sim_results[t_ind-1, bkt_col_index], sim_results[t_ind, bkt_col_index], waste_amounts)
                    

                            # if adsorber is full, remove from list of adsorbers
                            full_adsorber = adsorbers.pop(0)
                            # reset the states of the adsorber to zero (i.e., regenerate the adsorber)
                            full_adsorber.reset_states(t)  # Pass the current time to the reset_states method
                            # Place the adsorber in RR and determine completion time
                            full_adsorber.inRR = True
                            full_adsorber.set_doneRRtime(t+RRtime) 
                            # append the adsorber to the end of the list
                            adsorbers.append(full_adsorber)
                            switch_times.append(t)
                            col_switches += 1
                            # Update system configuration per the new adsorber order
                            valve1.dest = adsorbers[0]
                            adsorbers[0].source = valve1
                            adsorbers[0].dest = pump2
                            pump2.source = adsorbers[0]

                            valve2.dest = adsorbers[1]
                            adsorbers[1].source = valve2
                            adsorbers[1].dest = pump3
                            pump3.source = adsorbers[1]

                            valve3.dest = adsorbers[2]
                            adsorbers[2].source = valve3
                            adsorbers[2].dest = pump4
                            pump4.source = adsorbers[2]

                            valve4.dest = adsorbers[3]
                            adsorbers[3].source = valve4

                        # If not full yet, check for breakthrough and calculate waste                        
                        else:
                            # Use dictionary to retrieve column index in sim_results for the breakthrough receiving adsorber
                            if adsorbers[1].inRR == False:
                                # If it is available, breakthrough is received by adsorbers[1]
                                bkt_col_index = adsorber_to_column[adsorbers[1].name]
                                # Set the stage array to 1 for receiving breakthrough
                                Ads_stage_array[t_ind, bkt_col_index - 1] = 1
                                                        
                            elif adsorbers[1].inRR == True:
                                # If not available, the product exiting adsorbers[0] is waste
                                bkt_col_index = adsorber_to_column[adsorbers[0].name]

                            # If breakthough in adsorber, determine the area under the curve for that time interval (waste)
                            if sim_results[t_ind, bkt_col_index] >= 0.00001:
                                waste_calc(t, t+delta_t, sim_results[t_ind-1, bkt_col_index], sim_results[t_ind, bkt_col_index], waste_amounts)
                    
                    
                    else:
                        # If no loading:
                        waste_amounts.append(0)

                    
                    # Prepare initial conditions for the next step using updated adsorber states
                    conc_ads_bulk_bound_0 = create_initial_conditions_vector(adsorbers, K) 
                    # Reassign the initial conditions for Tank1 and Tank0 for next step
                    conc_vol_tank1_0 = [tank1concs[-1], tank1vols[-1]]             
                   
                    # C_capacity = C_capacity_func(100,col_switches,adsorbers)

                # Create a dotted line for the threshold concentration
                threshline = np.ones(len(t_array))
                threshline = threshline*conc*thresh
                                
                #%% Plot the BTC results
                inlet_conc = np.zeros(len(t_array))
                inlet_flow = np.zeros(len(t_array))

                for t_idx,t in enumerate(t_array):
                    inlet_conc[t_idx] = inlet_feed_concentration(t, conc, conc_amp)
                    inlet_flow[t_idx] = inlet_flow_function(t, flow, flow_amp)

                plt.figure(dpi=200)
                for i in range(len(adsorbers)):

                    plt.plot(sim_results[:,0],sim_results[:,i+1])
                plt.plot(sim_results[:,0],tank1concs[1:])
                plt.plot(sim_results[:,0],tank0concs[1:])
                plt.plot(sim_results[:,0], threshline, linestyle='dotted', color='grey')
                # plt.plot(sim_results[:,0],Ads_stage_array[:,i])
                plt.xlabel('Run Time (min)')
                plt.ylabel('Outlet Conc. (mg/mL)')
                plt.ylim([-2,4])         
                
                labels = ['Adsorber 1', 'Adsorber 2', 'Adsorber 3', 'Adsorber 4', 'Surge Vessel Conc.', 'Feed Conc.', 'Threshold Conc.', 'Harvest Flow Rate']
                plt.legend(labels, loc='lower left')

                plt.tight_layout()
                plt.show(block=False)


                # Also plot the volume of Tank0 and Tank1
                plt.figure(dpi=200)
                plt.plot(sim_results[:,0], tank1vols[1:])
                plt.plot(sim_results[:,0], tank0vols[1:])
                plt.title(f'Q = {flow}, C = {conc}, T = {thresh}')
                plt.xlabel('Run Time (min)')
                plt.ylabel('Tank Volume (mL)')
                # Create a secondary y-axis
                ax1 = plt.gca()
                ax2 = ax1.twinx()
                ax2.plot(sim_results[:, 0], inlet_flow, color='black', linestyle='dashed')
                ax2.set_ylabel('Flow Rate (mL/min)')
               # Set the same y-limits and y-ticks for all twin axes
                ax2.set_ylim([-0.2, 3.2])
                ax2.set_yticks([0, 1, 2, 3])
                lines = ax1.get_lines() + ax2.get_lines()
                labels = ['Surge Vessel Volume', 'Bioreactor Volume', 'Harvest Flow Rate']
                plt.legend(lines, labels, loc='upper left')
                plt.tight_layout()
                plt.show(block=False)


                #%% Analyze the switch time and histogram plot
                print("Flow: %.2f" % flow)
                print("Concentration: %.2f" % conc)
                print("Threshold: %.2f" % switch_thresh) 
                print("Regeneration Time: %i" % RRtime)
                                
                switch_num = len(switch_times)
                print("Number of switches (and resin volume used): %i" % switch_num)

                # Calculate total waste and total inlet mass from the integrated values*flow rate
                # Waste from breakthrough:
                wastes[flowind, concind, thresind] = sum(waste_amounts)*flow
                # Mass fed to adsorbers: 
                masses[flowind, concind, thresind] = sum(mass_amounts)*flow
                print('Mass fed to adsorbers: %.2f' % masses[flowind,concind,thresind])
                # Total mass into surge vessel:
                total_masses[flowind, concind, thresind] = sum(total_mass_amounts)*flow
                print('Total mass into surge vessel: %.2f' % total_masses[flowind, concind, thresind])
                
                # Determine total waste from overflow
                overflow_wastes[flowind, concind, thresind] = sum(overflow_waste_amounts)
                print('Overflow: %.2f mg' % overflow_wastes[flowind, concind, thresind])
                total_waste[flowind, concind, thresind] = wastes[flowind, concind, thresind] +  overflow_wastes[flowind, concind, thresind]
                print('Total waste: %.2f mg' % total_waste[flowind, concind, thresind])
                # Mass adsorbed = masses into adsorbers - waste by adsorbers
                mass_adsorbed = masses[flowind, concind, thresind] - wastes[flowind, concind, thresind]
                print('Mass adsorbed: %.2f' % mass_adsorbed)


                prods[flowind, concind, thresind] = mass_adsorbed / (num_cols*col_vols[flowind, concind]) / (time_steps/60) #mg/mL/hour
                print("Productivity at a threshold of %.2f: %.1f mg/mL/hr" % (switch_thresh,prods[flowind, concind, thresind]))
              
                # Waste percent = total waste by system/mass adsorbed by adsorbers
                waste_percent[flowind, concind, thresind] = total_waste[flowind, concind, thresind]*100 / mass_adsorbed
                print("Waste percentage: %.2f" % waste_percent[flowind, concind, thresind])
                percent_overflow_waste[flowind, concind, thresind] = overflow_wastes[flowind, concind, thresind]*100/total_waste[flowind, concind, thresind]
                print("Percent of waste from overflow: %.2f" % percent_overflow_waste[flowind, concind, thresind])
                
                
                total_runs = len(flows)*len(concs)*len(threshs)
                runs += 1
                print(tank1vols)
                print('Completed: %i of %i\n\n' % (runs, total_runs))
                

                # Plot the cycles
                fig, axs = plt.subplots(4,1, figsize=(10,9))
                for i in range(len(adsorbers)):

                    ax = axs[i]
    
                    # Plot each adsorber's concentration profile and stage sequencing
                    ax.plot(sim_results[:, 0], sim_results[:, i+1], label='Conc. (mg/mL)')
                    ax.set_xlabel('Run Time (min)')
                    ax.set_ylabel('Conc. (mg/mL)')
                    
                    # Create a twin axis
                    twin_ax = ax.twinx()
                    twin_ax.plot(sim_results[:, 0], Ads_stage_array[:, i], color='orange', label='Adsorber Stage')
                    twin_ax.set_ylabel('Adsorber Stage')

                    # Set the same y-limits and y-ticks for all twin axes
                    twin_ax.set_ylim([-1.2, 2.2])
                    twin_ax.set_yticks([-1, 0, 1, 2])
                    
                fig.tight_layout()
                plt.show()

                # Raise expection for when the surge vessel dries up or overflows
                for t_ind, time in enumerate(t_array):
                    if tank1vols[t_ind] < 0:
                        print(f'\n\nSURGE VESSEL DRIED UP (V <= =) ... time of error: {t_array[t_ind-1]}\n\n')

                    if tank1vols[t_ind] > tank1.max_volume:
                        print(f'\n\nSURGE VESSEL OVERFLOWED (V > max) ... time of error: {t_array[t_ind-1]}\n\n')


                # # Calculate the differences between consecutive switch times
                # differences = [switch_times[i+1] - switch_times[i] for i in range(len(switch_times)-1)]

                # # Convert the differences to a numpy array
                # differences_array = np.array(differences)

                # # Function to plot histogram with a given bin width
                # def plot_histogram(data, bin_width):
                #     # Define the bins
                #     min_bin = min(data)
                #     max_bin = max(data) + bin_width  # Add bin_width to include the max value in the histogram
                #     bins = np.arange(min_bin, max_bin, bin_width)

                #     # Create a histogram for the data
                #     plt.hist(data, bins=bins, edgecolor='black')

                #     # Set the ticks on the x-axis to match the bins
                #     plt.xticks(bins)

                #     # Ensure y-axis has only integer ticks
                #     plt.gca().yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

                #     # Add title and labels as needed
                #     plt.title('Histogram of Differences in Switch Times, %.2f Threshold' % switch_thresh)
                #     plt.xlabel('Difference Values')
                #     plt.ylabel('Frequency')

                #     # Show the plot
                #     plt.show()

                # # Example usage
                # bin_width = 2  # You can change this value to adjust the bin width
                # plot_histogram(differences_array, bin_width)

 
# export the waste and productivity values to Excel to organize
if len(threshs) == 1:
    flat_waste = waste_percent.reshape(len(flows), len(concs))
    flat_prods = prods.reshape(len(flows), len(concs))
elif len(concs) == 1: 
    flat_waste = waste_percent.reshape(len(flows), len(threshs))
    flat_prods = prods.reshape(len(flows), len(threshs))
elif len(flows) == 1: 
    flat_waste = waste_percent.reshape(len(concs), len(threshs))
    flat_prods = prods.reshape(len(concs), len(threshs))

# Convert the product and waste results to data frames to export
prods_data = pd.DataFrame(flat_prods)
wastes_data = pd.DataFrame(flat_waste)
volume_data = pd.DataFrame(col_vols)

# The Excel files must be closed to overwrite the data from the previous trial
wastes_data.to_excel('Regen Time Tank 4-Col Run Wastes.xlsx', index = False, header = False)
prods_data.to_excel('Regen Time Tank 4-Col Run Prods.xlsx', index = False, header = False)
volume_data.to_excel('Regen Time Tank 4-Column Volumes.xlsx', index = False, header = False)
