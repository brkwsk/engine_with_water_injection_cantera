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

#space for storage of the input variables
#ca_values = np.array([])
ca_values = []
#wf_values = np.array([])
wf_values = []
input_values = {"ca":ca_values,"wf":wf_values}


class GUI_Input:
    """This class provides a GUI, takes input valuesof crank angle at the beginning 
    of water injection and water/fuel ratio, and returns them as a dictionary containing arrays. """
    
    def __init__(self,master):
        
        self.master = master
        
        #make frames for content of GUI
        frame1 = tk.LabelFrame(master, text = "Current values",relief = "groove", borderwidth = 5)
        frame1.pack(fill = tk.BOTH,side ="left", expand = True)
        
        frame2 = tk.LabelFrame(master,text = "Modify variables",relief = "groove", borderwidth = 5)
        frame2.pack(fill = tk.BOTH,side ="right", expand = True)
        
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
        lambda: self.add_to_entry(wf_values, self.current_wf, self.new_wf_entry.get()))
        self.wf_add_button.pack(fill = tk.BOTH, side = "left", expand = True)
        
        self.wf_remove_button = tk.Button(frame22, text = "Remove", command = 
        lambda: self.remove_from_entry(wf_values, self.current_wf, self.new_wf_entry.get()))
        self.wf_remove_button.pack(fill = tk.BOTH, side = "right", expand = True)
        
        self.new_ca_label = tk.Label(frame23, text = "CA")
        self.new_ca_label.pack(fill = tk.BOTH, side = "left", expand = True)
        
        self.new_ca_entry = tk.Entry(frame23)
        self.new_ca_entry.pack(fill = tk.BOTH, side = "right", expand = True)
        
        #command functions for the buttons are defined below __init__
        self.ca_add_button = tk.Button(frame24, text = "Add", command = 
        lambda: self.add_to_entry(ca_values, self.current_ca, self.new_ca_entry.get()))
        self.ca_add_button.pack(fill = tk.BOTH, side = "left", expand = True)
        
        self.ca_remove_button = tk.Button(frame24, text = "Remove", command = 
        lambda: self.remove_from_entry(ca_values, self.current_ca, self.new_ca_entry.get()))
        self.ca_remove_button.pack(fill = tk.BOTH, side = "right", expand = True)
        
        self.run_button = tk.Button(frame2, text = "Run simulations")
        self.run_button.pack(fill= tk.BOTH, side = "bottom")
        
        #functions for GUI buttons
    def add_to_entry(self, container, entry, value):
        #container = np.append(container, value)
        container.append(value)
        entry.delete(0, tk.END)
        entry.insert(0, container)
        #print(container)
        
    def remove_from_entry(self, container, entry, value):
        #container = np.delete(container, np.where(container == value))
        container.remove(value)
        entry.delete(0, tk.END)
        entry.insert(0, container)
        #print(container)

        
       
    
#initialize GUI   
input_window = tk.Tk()
app = GUI_Input(input_window)
input_window.wm_title("Values of variables for smulations")
#show GUI
input_window.mainloop()

        

