# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 17:57:50 2021

@author: Andrzej

This program takes input values of crank angle (CA) 
and water/fuel ratio (w/f) through GUI,
runs simulation for each combination of input values,
and saves output as text files containing parameters of
every 10th step of simulation plus p-V plots. 
"""
#import cantera
import tkinter as tk
#from matplotlib import pyplot as plt
#import numpy as np

class GUI_Input:
    """This class provides a GUI, takes input values of crank angle at the beginning 
    of water injection and water/fuel ratio, and saves them in "input_values" dictionary 
    containing lists of values."""
    
    def __init__(self,master):
        
        self.master = master
        
        #space for storage of the input variables
        self.input_values = {"ca":[],"wf":[]} 
        
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
        
        #make entry fields, buttons and labels
        self.current_wf_label = tk.Label(frame1, text = "Water/fuel ratio", relief = "ridge", borderwidth = 2)
        self.current_wf_label.pack(fill = tk.BOTH, expand = True)
        
        self.current_wf = tk.Entry(frame1)
        self.current_wf.pack(fill = tk.BOTH, expand = True)
        
        self.current_ca_label = tk.Label(frame1,text = "Crank angle at the beginning\nof water injection", relief = "ridge", borderwidth = 2)
        self.current_ca_label.pack(fill = tk.BOTH, expand = True)
        
        self.current_ca = tk.Entry(frame1)
        self.current_ca.pack(fill = tk.BOTH,expand = True)
        
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
                self.messages.insert(tk.INSERT,"Running simulation: CA = {}, w/f = {}\n".format(ca,wf))
                Simulate(wf,ca)



class Simulate:
    """This class manages simulation and output for a single set of input variables.
    At this moment it is a placeholder"""
    
    def __init__(self,wf,ca):
        self.wf = wf
        self.ca = ca
    pass
        



#space for storage of the input variables
#ca_values = np.array([])
ca_values = []
#wf_values = np.array([])
wf_values = []
input_values = {"ca":ca_values,"wf":wf_values}        
       
    
#initialize GUI   
input_window = tk.Tk()
app = GUI_Input(input_window)
input_window.wm_title("Values of variables for smulations")
#show GUI
input_window.geometry("600x180")
input_window.mainloop()

        

