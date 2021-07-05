# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 17:57:50 2021

@author: Andrzej

This program takes input values of crank angle (CA) 
and water/fuel ratio (w/f) through GUI,
runs simulation for each combination of input values,
and saves output as text files containing parameters of
every 10th step of simulation plus p-V plots.
It is based on Cantera example ic_engine.py. 
"""
import cantera as ct
import tkinter as tk
from matplotlib import pyplot as plt
import numpy as np
from scipy.integrate import trapz

class GUI_Input:
    """This class provides a GUI, takes input values of crank angle at the beginning 
    of water injection, pressure of injection and water/fuel ratio, and saves them in "input_values" dictionary 
    containing lists of values."""
    
    def __init__(self,master):
        
        self.master = master
        
        #space for storage of the input variables
        self.input_values = {"ca":[],"wf":[],"p":[]} 
        
        #make frames for content of GUI
        #frame1 (left) shows the submitted values
        frame1 = tk.LabelFrame(master, text = "Current values",relief = "groove", borderwidth = 5)
        frame1.pack(fill = tk.BOTH, side = "left", expand = True)
        
        #frame2 (middle) manages editing values
        frame2 = tk.LabelFrame(master,text = "Modify variables",relief = "groove", borderwidth = 5)
        frame2.pack(fill = tk.BOTH, side = "left", expand = True)
        
        #frame3 (right) contains a message window
        frame3 = tk.LabelFrame(master, text = "Message window", relief = "groove", borderwidth = 5)
        frame3.pack(fill = tk.BOTH, side = "left", expand = True)
        
        #sub-frames allow vertical and horizontal management of widgets
        frame21 = tk.LabelFrame(frame2)
        frame21.pack(fill = tk.BOTH,expand = True)
        frame22 = tk.LabelFrame(frame2)
        frame22.pack(fill = tk.BOTH,expand = True)
        frame23 = tk.LabelFrame(frame2)
        frame23.pack(fill = tk.BOTH,expand = True)
        frame24 = tk.LabelFrame(frame2)
        frame24.pack(fill = tk.BOTH,expand = True)
        frame25 = tk.LabelFrame(frame2)
        frame25.pack(fill = tk.BOTH,expand = True)
        frame26 = tk.LabelFrame(frame2)
        frame26.pack(fill = tk.BOTH,expand = True)
        
        #make entry fields, buttons and labels
        self.current_wf_label = tk.Label(frame1, text = "Water/fuel ratio", relief = "ridge", borderwidth = 2)
        self.current_wf_label.pack(fill = tk.BOTH, expand = True)
        
        self.current_wf = tk.Entry(frame1)
        self.current_wf.pack(fill = tk.BOTH, expand = True)
        
        self.current_ca_label = tk.Label(frame1,text = "Crank angle at the beginning\nof water injection [degrees]", relief = "ridge", borderwidth = 2)
        self.current_ca_label.pack(fill = tk.BOTH, expand = True)
        
        self.current_ca = tk.Entry(frame1)
        self.current_ca.pack(fill = tk.BOTH,expand = True)
        
        self.current_pressure_label = tk.Label(frame1,text = "Pressure of injection [bar]", relief = "ridge", borderwidth = 2)
        self.current_pressure_label.pack(fill = tk.BOTH, expand = True)
        
        self.current_pressure = tk.Entry(frame1)
        self.current_pressure.pack(fill = tk.BOTH, expand = True)
        
        self.new_wf_label = tk.Label(frame21, text = "w/f")
        self.new_wf_label.pack(fill = tk.BOTH, side = "left", expand = True)
        
        self.new_wf_entry = tk.Entry(frame21)
        self.new_wf_entry.pack(fill = tk.BOTH, side = "right", expand = True)
        
        #command functions for the buttons are defined below __init__
        self.wf_add_button = tk.Button(frame22, text = "Add", command = 
        lambda: self.add_to_entry(self.input_values["wf"], self.current_wf, self.new_wf_entry))
        self.wf_add_button.pack(fill = tk.BOTH, side = "left", expand = True)
        
        self.wf_remove_button = tk.Button(frame22, text = "Remove", command = 
        lambda: self.remove_from_entry(self.input_values["wf"], self.current_wf, self.new_wf_entry))
        self.wf_remove_button.pack(fill = tk.BOTH, side = "right", expand = True)
        
        self.new_ca_label = tk.Label(frame23, text = "CA")
        self.new_ca_label.pack(fill = tk.BOTH, side = "left", expand = True)
        
        self.new_ca_entry = tk.Entry(frame23)
        self.new_ca_entry.pack(fill = tk.BOTH, side = "right", expand = True)
        
        #command functions for the buttons are defined below __init__
        self.ca_add_button = tk.Button(frame24, text = "Add", command = 
        lambda: self.add_to_entry(self.input_values["ca"], self.current_ca, self.new_ca_entry))
        self.ca_add_button.pack(fill = tk.BOTH, side = "left", expand = True)
        
        self.ca_remove_button = tk.Button(frame24, text = "Remove", command = 
        lambda: self.remove_from_entry(self.input_values["ca"], self.current_ca, self.new_ca_entry))
        self.ca_remove_button.pack(fill = tk.BOTH, side = "right", expand = True)
        
        self.new_pressure_label = tk.Label(frame25, text = "Pressure")
        self.new_pressure_label.pack(fill = tk.BOTH, side = "left", expand = True)
        
        self.new_pressure_entry = tk.Entry(frame25)
        self.new_pressure_entry.pack(fill = tk.BOTH, side = "right", expand = True)
        
        #command functions for the buttons are defined below __init__
        self.pressure_add_button = tk.Button(frame26, text = "Add", command = 
        lambda: self.add_to_entry(self.input_values["p"], self.current_pressure, self.new_pressure_entry))
        self.pressure_add_button.pack(fill = tk.BOTH, side = "left", expand = True)
        
        self.pressure_remove_button = tk.Button(frame26, text = "Remove", command = 
        lambda: self.remove_from_entry(self.input_values["p"], self.current_pressure, self.new_pressure_entry))
        self.pressure_remove_button.pack(fill = tk.BOTH, side = "right", expand = True)
        
        self.run_button = tk.Button(frame2, text = "Run simulations", command = 
        lambda: self.run_simulations())
        self.run_button.pack(fill= tk.BOTH, side = "bottom")
        
        self.messages = tk.Text(frame3)
        self.messages.pack(fill = tk.BOTH, expand = True)
        
    #functions for GUI buttons
    def add_to_entry(self, container, entry_out, entry_in):
        if entry_in.get() == "":        #check if there is an input
            self.messages.insert(tk.INSERT, "No value entered" + "\n")
        else:
            value = float(entry_in.get())
            if value not in container:      #check if the value has already been entered
                container.append(value)     #update the container
                container.sort()
                entry_out.delete(0, tk.END)     #update output entry
                entry_out.insert(0, container)                
            else:
                self.messages.insert(tk.INSERT, "Value already added" + "\n")
            entry_in.delete(0, tk.END)      #clear input entry
        
    def remove_from_entry(self, container, entry_out, entry_in):
        if entry_in.get() == "":        #check if there is an input
            self.messages.insert(tk.INSERT, "No value entered" + "\n")
        else:
            value = float(entry_in.get())
            if value in container:      #check if the value has been entered
                container.remove(value)     #update container
                entry_out.delete(0, tk.END)     #update output entry
                entry_out.insert(0, container)              
            else:
                self.messages.insert(tk.INSERT, "Value not in set" + "\n")
            entry_in.delete(0, tk.END)      #clear input entry

    #running simulations while iterating through variable values
    def run_simulations(self):
        for wf in self.input_values["wf"]:
            for ca in self.input_values["ca"]:
                for p in self.input_values["p"]:
                    self.messages.insert(tk.INSERT,"Running simulation: CA = {} degrees, w/f = {}, p = {} bar\n".format(ca,wf,p))
                    Simulate(wf,ca,p)



class Simulate:
    """This class manages simulation and output for a single set of input variables. It is a modified
    version of Cantera example ic_engine.py, available at Cantera website.
    At this moment it is a placeholder"""
    
    def __init__(self,wf,ca,p):
        self.wf = wf
        self.deg_water_open = ca
        self.p_water = p*1e5 # Pa
        
        #########################################################################
        # Input Parameters
        #########################################################################

        # reaction mechanism, kinetics type and compositions
        reaction_mechanism = 'blanquart.cti'
        water_injection_mechanism = "liquidvapor.cti" # required to include evaporation of water
        phase_name = 'I-C8H18'
        comp_air = 'o2:1, n2:3.76'
        comp_fuel = 'I-C8H18:1' # approximation of gasoline composition
        comp_spark = 'OH:1' # highly reactive radicals allow reduction of the mass of spark
        
        # engine parameters
        self.f = 3000. / 60.  # engine speed [1/s] (3000 rpm)
        self.stroke = 0.15 # stroke [m]
        
        epsilon = 15.  # compression ratio [-]
        d_piston = 0.06  # piston diameter [m]
        
        A_piston = np.pi * d_piston ** 2 / 4 # area of piston [m**2]
        V_H = self.stroke * A_piston  # displaced volume [m**3]
        V_oT = V_H / (epsilon - 1.) # minimal volume [m**3]
        
        # parameters of spark ignition
        e_spark = 50. # approximateenergy of spark [J]
        T_spark = 10000. # temperature of spark [K]
        
        # turbocharger temperature, pressure, and composition
        T_inlet = 300.  # K
        p_inlet = 1.4e5  # Pa
        comp_inlet = comp_air

        # outlet pressure
        p_outlet = 1.0e5  # Pa

        # fuel properties (gaseous!)
        T_injector = 300.  # K
        p_injector = 200.0e5  # Pa
        comp_injector = comp_fuel
        
        # water properties
        T_water = 300. # K

        # ambient properties
        T_ambient = 300.  # K
        p_ambient = 1e5  # Pa
        comp_ambient = comp_air

        # Inlet valve friction coefficient, open and close timings
        inlet_valve_coeff = 1.e-6
        inlet_open = self.rad(-18.)
        inlet_close = self.rad(198.)

        # Outlet valve friction coefficient, open and close timings
        outlet_valve_coeff = 1.e-6
        outlet_open = self.rad(522.)
        outlet_close = self.rad(18.)

        # Fuel mass, injector open and close timings
        injector_open = self.rad(345.)
        injector_close = self.rad(360.)
        lmbd = 1.2 # air-fuel ratio
        V_air = 0.44927e-3 # volume of air at the moment of inlet valve closing, measured in simulation [m**3]
        p_air = 1.1e5 # in-cyllinder pressure at the moment of inlet valve closing, measured in simulation [m**3]
        T_air = 330. # in_ cyllinder temperature at the moment of inlet valve closing, measured in simulation [m**3]
        stoichiometry = 12.5 # stoichiometric mole proportion between oxygen and isooctane
        o2_in_air = 4.76
        n_fuel = p_air*V_air/lmbd/12.5/stoichiometry/8.314/T_air/o2_in_air # moles of fuel
        injector_mass = n_fuel * 114./1000  # kg
        
        # Ignition timing
        start_spark = self.rad(355.)
        end_spark = self.rad(356.)
        
        # Water mass and injector timings
        # Injection is prolonged to simulate evaporation rate
        self.water_open = self.rad(self.deg_water_open)
        self.water_mass = injector_mass * self.wf
        self.water_evaporation_time = self.calculate_water_evaporation_time() # s
        self.water_close = self.water_open + self.crank_angle(self.water_evaporation_time)
        
        second_water_evaporation_time = 1.
        self.water_mass_2 = 0.
        water_close_2 = start_spark
        if self.water_close > start_spark: # water has not evaporated before ignition - increase evaporation rate due to higher temperatures
            self.water_close = start_spark
            second_water_evaporation_time = self.calculate_second_water_evaporation_time()
            water_close_2 = self.water_close + self.crank_angle(second_water_evaporation_time)
            print("Water evaporated after ignition")
            if water_close_2 > outlet_open: # if water didn't evaporate completely: artificially shorten evaporation
                print("Water didn't fully evaporate")
                second_water_evaporation_time = self.time(outlet_open - self.water_close)
                self.water_mass_2 = second_water_evaporation_time * self.evaporation_rate_2/1000.
                water_close_2 = outlet_open
        
        # Simulation time and parameters
        sim_n_revolutions = 6   # 3 cycles of engine work
        delta_T_max = 20.
        rtol = 1.e-12
        atol = 1.e-16
        
        #####################################################################
        # Set up Reactor Network
        #####################################################################
        
        # load reaction mechanism
        gas = ct.Solution(reaction_mechanism, id = phase_name)
        w = ct.PureFluid(water_injection_mechanism,"water") # different mechanism for water injector allows calculations using evaporation heat
        
        # define initial state and set up reactor
        gas.TPX = T_inlet, p_inlet, comp_inlet
        cyllinder = ct.IdealGasReactor(gas)
        cyllinder.volume = V_oT
        
        # define inlet state
        inlet = ct.Reservoir(gas)
        
        # inlet valve
        inlet_valve = ct.Valve(inlet, cyllinder)
        inlet_delta = np.mod(inlet_close - inlet_open, 4 * np.pi)
        inlet_valve.valve_coeff = inlet_valve_coeff
        inlet_valve.set_time_function(
            lambda t: np.mod(self.crank_angle(t) - inlet_open, 4 * np.pi) < inlet_delta)
        
        # define fuel injector state (gaseous!)
        gas.TPX = T_injector, p_injector, comp_injector
        fuel_injector = ct.Reservoir(gas)
        
        # define water injector state
        w.TP = T_water, 1.0e5
        water_injector = ct.Reservoir(w)
        
        # define spark state
        gas.TPX = T_spark, p_injector, comp_spark
        igniter = ct.Reservoir(gas)
        
        # fuel injector is modeled as a mass flow controller
        fuel_injector_mfc = ct.MassFlowController(fuel_injector, cyllinder)
        fuel_injector_delta = np.mod(injector_close - injector_open, 4 * np.pi)
        fuel_injector_t_open = self.time(injector_close - injector_open)
        fuel_injector_mfc.mass_flow_coeff = injector_mass / fuel_injector_t_open
        fuel_injector_mfc.set_time_function(
            lambda t: np.mod(self.crank_angle(t) - injector_open, 4 * np.pi) < fuel_injector_delta)
        
        # water injector is modeled as two mass flow controllers, activating before and after ignition
        water_injector_mfc = ct.MassFlowController(water_injector, cyllinder)
        water_injector_mfc_2 = ct.MassFlowController(water_injector, cyllinder)
        
        water_injector_delta = np.mod(self.water_close - self.water_open, 4 * np.pi)
        water_injector_delta_2 = np.mod(water_close_2 - self.water_close, 4 * np.pi)
        
        water_injector_mfc.mass_flow_coeff = self.water_mass / self.water_evaporation_time
        water_injector_mfc.set_time_function(
            lambda t: np.mod(self.crank_angle(t) - self.water_open, 4 * np.pi) < water_injector_delta)
        
        water_injector_mfc_2.mass_flow_coeff = self.water_mass_2 / second_water_evaporation_time
        water_injector_mfc_2.set_time_function(
            lambda t: np.mod(self.crank_angle(t) - self.water_close, 4 * np.pi) < water_injector_delta_2)
        
        # spark ignition is modeled as a mass flow controller (approximation!)
        spark_mfc = ct.MassFlowController(igniter, cyllinder)
        spark_delta = np.mod(end_spark - start_spark, 4 * np.pi)
        spark_t_open = self.time(end_spark - start_spark)
        spark_mass = e_spark/igniter.thermo.cp_mass/T_spark
        spark_mfc.mass_flow_coeff = spark_mass / spark_t_open
        spark_mfc.set_time_function(
            lambda t: np.mod(self.crank_angle(t) - start_spark, 4 * np.pi) < spark_delta)
        
        # define outlet pressure (temperature and composition don't matter)
        gas.TPX = T_ambient, p_outlet, comp_ambient
        outlet = ct.Reservoir(gas)
        
        # outlet valve
        outlet_valve = ct.Valve(cyllinder, outlet)
        outlet_delta = np.mod(outlet_close - outlet_open, 4 * np.pi)
        outlet_valve.valve_coeff = outlet_valve_coeff
        outlet_valve.set_time_function(
            lambda t: np.mod(self.crank_angle(t) - outlet_open, 4 * np.pi) < outlet_delta)
        
        # define ambient pressure (temperature and composition don't matter)
        gas.TPX = T_ambient, p_ambient, comp_ambient
        ambient_air = ct.Reservoir(gas)
        
        # piston is modeled as a moving wall
        piston = ct.Wall(ambient_air, cyllinder)
        piston.area = A_piston
        piston.set_velocity(self.piston_speed)
        
        # create a reactor network containing the cylinder and limit advance step
        sim = ct.ReactorNet([cyllinder])
        sim.rtol, sim.atol = rtol, atol
        cyllinder.set_advance_limit('temperature', delta_T_max)
        
        #####################################################################
        # Run Simulation
        #####################################################################
        
        # set up output data arrays
        states = ct.SolutionArray(
            cyllinder.thermo,
            extra=('t', 'ca', 'V', 'm', 'mdot_in', 'mdot_out', 'dWv_dt'),
        )
        
        # simulate with a maximum resolution of 1 deg crank angle
        dt = 1. / (360 * self.f)
        t_stop = sim_n_revolutions / self.f
        while sim.time < t_stop:
        
            # perform time integration
            sim.advance(sim.time + dt)
        
            # calculate results to be stored
            dWv_dt = - (cyllinder.thermo.P - ambient_air.thermo.P) * A_piston * \
                self.piston_speed(sim.time)
        
            # append output data
            states.append(cyllinder.thermo.state,
                          t=sim.time, ca=self.crank_angle(sim.time),
                          V=cyllinder.volume, m=cyllinder.mass,
                          mdot_in=inlet_valve.mass_flow_rate,
                          mdot_out=outlet_valve.mass_flow_rate,
                          dWv_dt=dWv_dt)
        
        #######################################################################
        # Plot Results in matplotlib
        #######################################################################
        
        def ca_ticks(t):
            """Helper function converts time to rounded crank angle."""
            return np.round(self.crank_angle(t) * 180 / np.pi, decimals=1)
        
        t = states.t
        
        # pressure and temperature
        xticks = np.arange(0, 0.18, 0.02)
        fig, ax = plt.subplots(nrows=2)
        ax[0].plot(t, states.P / 1.e5)
        ax[0].set_ylabel('$p$ [bar]')
        ax[0].set_xlabel(r'$\phi$ [deg]')
        ax[0].set_xticklabels([])
        ax[1].plot(t, states.T)
        ax[1].set_ylabel('$T$ [K]')
        ax[1].set_xlabel(r'$\phi$ [deg]')
        ax[1].set_xticks(xticks)
        ax[1].set_xticklabels(ca_ticks(xticks))
        plt.show()
        
        # p-V diagram
        fig, ax = plt.subplots()
        ax.plot(states.V[t > 0.04] * 1000, states.P[t > 0.04] / 1.e5)
        ax.set_xlabel('$V$ [l]')
        ax.set_ylabel('$p$ [bar]')
        plt.show()
        
        # T-S diagram
        fig, ax = plt.subplots()
        ax.plot(states.m[t > 0.04] * states.s[t > 0.04], states.T[t > 0.04])
        ax.set_xlabel('$S$ [J/K]')
        ax.set_ylabel('$T$ [K]')
        plt.show()
        
        # heat of reaction and expansion work
        fig, ax = plt.subplots()
        ax.plot(t, 1.e-3 * states.heat_release_rate * states.V, label=r'$\dot{Q}$')
        ax.plot(t, 1.e-3 * states.dWv_dt, label=r'$\dot{W}_v$')
        ax.set_ylim(-1e2, 1e3)
        ax.legend(loc=0)
        ax.set_ylabel('[kW]')
        ax.set_xlabel(r'$\phi$ [deg]')
        ax.set_xticks(xticks)
        ax.set_xticklabels(ca_ticks(xticks))
        plt.show()
        
        # gas composition
        fig, ax = plt.subplots()
        ax.plot(t, states('o2').X, label='O2')
        ax.plot(t, states('co2').X, label='CO2')
        ax.plot(t, states('co').X, label='CO')
        ax.plot(t, states('i-c8h18').X * 10, label='iso-octane x10')
        ax.legend(loc=0)
        ax.set_ylabel('$X_i$ [-]')
        ax.set_xlabel(r'$\phi$ [deg]')
        ax.set_xticks(xticks)
        ax.set_xticklabels(ca_ticks(xticks))
        plt.show()
        
        ######################################################################
        # Integral Results
        ######################################################################
        
        # heat release
        Q = trapz(states.heat_release_rate * states.V, t)
        output_str = '{:45s}{:>4.1f} {}'
        print(output_str.format('Heat release rate per cylinder (estimate):',
                                Q / t[-1] / 1000., 'kW'))
        
        # expansion power
        W = trapz(states.dWv_dt, t)
        print(output_str.format('Expansion power per cylinder (estimate):',
                                W / t[-1] / 1000., 'kW'))
        
        # efficiency
        eta = W / Q
        print(output_str.format('Efficiency (estimate):', eta * 100., '%'))
        
        # CO emissions
        MW = states.mean_molecular_weight
        CO_emission = trapz(MW * states.mdot_out * states('CO').X[:, 0], t)
        CO_emission /= trapz(MW * states.mdot_out, t)
        print(output_str.format('CO emission (estimate):', CO_emission * 1.e6, 'ppm'))
        
        
        #####################################################################
        # Set up IC engine Functions
        #####################################################################
        

    def crank_angle(self,t):
        """Convert time to crank angle"""
        return np.remainder(2 * np.pi * self.f * t, 4 * np.pi)
    
    def time(self, ca):
        """Convert crank angle to time"""
        return ca/ 2/ np.pi/ self.f
    
    def piston_speed(self,t):
        """Approximate piston speed with sinusoidal velocity profile"""
        return - self.stroke / 2 * 2 * np.pi * self.f * np.sin(self.crank_angle(t))
        pass
    
    def rad(self, degrees):
        """Convert degrees to radians"""
        return degrees * np.pi/180.
    
    def calculate_water_evaporation_time(self):
        T_evaporation = ((840. - (355. - self.deg_water_open)/2. *520./175.)+460.)/2. # weighed average temperature during evaporation of water
        self.evaporation_rate_650 = 0.0183 * self.p_water/1.0e5 + 0.375 # evaporation rate at 650 K [g/s]
        self.evaporation_rate = self.evaporation_rate_650 * (T_evaporation -300.)/(650.-300.) # g/s
        return self.water_mass * 1000./self.evaporation_rate
    
    def calculate_second_water_evaporation_time(self):
        T_evaporation_2 = (3130.+1550.)/2.
        self.evaporation_rate_2 = self.evaporation_rate_650 * (T_evaporation_2 - 300.)/(650.-300.)
        self.water_evaporation_time = self.time(self.water_close - self.water_open)
        self.water_mass_2 = self.water_mass - self.water_evaporation_time * self.evaporation_rate/1000.
        self.water_mass = self.water_mass - self.water_mass_2
        return self.water_mass_2 * 1000./self.evaporation_rate_2
       
    
#initialize GUI   
input_window = tk.Tk()
app = GUI_Input(input_window)
input_window.wm_title("Values of variables for smulations")
#show GUI
input_window.geometry("600x210")
input_window.mainloop()

        

