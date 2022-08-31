# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 19:27:32 2022

@author: Adrian
"""

from math import pi, exp
import numpy as np

## convergence: Function that tests the convergence considering successive
## iterations of the tank.
## Parameters:
##      - dict1: dictionary with temperatures of the previous iteration
##      - dict1: dictionary with temperatures of the current iteration
##      - tolerance: number to define how small the difference between the iterations can be
def convergence(dict1,dict2,tolerance):
    assert len(dict1) == len(dict2)
    N_nodes = len(dict1)
    List1 = [ abs((dict1['node_'+str(i)]['T_ave'] -
                   dict2['node_'+str(i)]['T_ave'])/(dict1['node_'+str(i)]['T_ave']+273.15)) < tolerance 
             for i in range(1, N_nodes+1) ]
    List2 = [ abs((dict1['node_'+str(i)]['T_final'] -
                   dict2['node_'+str(i)]['T_final'])/(dict1['node_'+str(i)]['T_final']+273.15)) < tolerance 
             for i in range(1, N_nodes+1) ]
    return all(List1) and all(List2)

## temperatures_ordered: function that receives the nodal temperatures of the tank and 
def temperatures_ordered(dict1):
    List1 = [ dict1[i] >= dict1[i+1] for i in range(1, len(dict1)) ]
    return all(List1)

def mix_nodes(dict1):
    mixed_nodes = [ {'original_node_numbers': [i], 'temperature': dict1[i], 'mass':1} for i in dict1 ]
    while True:
        List1 = [ mixed_nodes[i]['temperature'] >= mixed_nodes[i+1]['temperature'] for i in range(len(mixed_nodes) - 1) ]
        if all(List1):
            break
        for i in range(len(mixed_nodes) - 1):
            if mixed_nodes[i]['temperature'] < mixed_nodes[i+1]['temperature']:
                total_mass = mixed_nodes[i]['mass'] + mixed_nodes[i+1]['mass']
                temperature = ((mixed_nodes[i]['mass']/total_mass)*mixed_nodes[i]['temperature'] +
                               (mixed_nodes[i+1]['mass']/total_mass)*mixed_nodes[i+1]['temperature'])
                nodes_list = mixed_nodes[i]['original_node_numbers'] + mixed_nodes[i+1]['original_node_numbers']
                mixed_nodes[i] = {'original_node_numbers': nodes_list, 'temperature': temperature, 'mass': total_mass}
                mixed_nodes.pop(i+1)
                break
    R1 = {}
    for i in range(len(mixed_nodes)):
        for node_number in mixed_nodes[i]['original_node_numbers']:
            R1[node_number] = mixed_nodes[i]['temperature']
    return R1
    
class Storage_Tank:
    def __init__(self, height, volume, N_nodes, initial_temperature, top_loss_coeff,
                 edge_loss_coeff, bottom_loss_coeff, inlet_1_node, outlet_1_node,
                 inlet_2_node, outlet_2_node, fluid_density, fluid_spec_heat, fluid_thermal_conductivity):
        self.height = height ## height of the tank, m
        self.volume = volume ## volume of the tank, m^3
        self.N_nodes = N_nodes ## number of nodes
        self.T = {'node_'+str(i+1): initial_temperature for i in range(N_nodes)} ## set initial temperatures of the nodes
        self.top_loss_coeff = top_loss_coeff ## loss coefficient at the top
        self.edge_loss_coeff = edge_loss_coeff ## loss coefficient at the edge
        self.bottom_loss_coeff = bottom_loss_coeff ## loss coefficient at the bottom
        self.cross_Area = volume/height ## cross-sectional area
        self.radius = (self.cross_Area/pi)**0.5 ## radius of the tank 
        self.node_height = height/N_nodes ## height of each node
        self.nodal_edge_area = 2*pi*self.radius*self.node_height ## lateral area of each node
        self.node_volume = self.volume/N_nodes ## volume of each node
        inlets_and_outlets = {} ## dictionary to register the nodes where the inlets and outlets are
        for i in range(1,N_nodes+1):
            node_i_list = []
            if inlet_1_node == i:
                node_i_list.append('inlet_1')
            if outlet_1_node == i:
                node_i_list.append('outlet_1')
            if inlet_2_node == i:
                node_i_list.append('inlet_2')
            if outlet_2_node == i:
                node_i_list.append('outlet_2')
            inlets_and_outlets['node_'+str(i)] = node_i_list
        self.inlets_and_outlets = inlets_and_outlets
        self.fluid_density = fluid_density ## density of the fluid contained in the tank
        self.fluid_spec_heat = fluid_spec_heat ## specific heat of the fluid contained in the tank
        self.fluid_thermal_conductivity = fluid_thermal_conductivity ## thermal conductivity of the fluid
        self.node_thermal_capacity = self.node_volume*self.fluid_density*self.fluid_spec_heat ## thermal capacity of each node
    ## Tank.iteration: it executes a single iteration to determine the temperatures of the nodes at the end of
    ## the time step and the average temperatures during the time step
    ## It usually must be called more than one time to converge to definitve outputs given a certain input.
    def iteration(self, inlet_1_flow_rate, inlet_1_T, inlet_2_flow_rate, inlet_2_T,
                  T_out_top, T_out_edge, T_out_bottom, previous_result, time_step):
        Sol = {}
        for i in range(1, self.N_nodes + 1):
            flow_in = 0 ## flow coming from inlets to the node
            flow_out = 0 ## flow going out of the node through outlets of the tank
            ## values of a and b are computed, according to the mathematical reference of type 158 of TRNSYS
            a = -self.nodal_edge_area*self.edge_loss_coeff
            b = self.nodal_edge_area*self.edge_loss_coeff*T_out_edge
            ## if the node has inlets or outlets, flow_in and flow_out are updated
            ## this also modifies the values of a and b
            if 'inlet_1' in self.inlets_and_outlets['node_'+str(i)]:
                flow_in  = flow_in + inlet_1_flow_rate
                a = a - inlet_1_flow_rate*self.fluid_spec_heat
                b = b + inlet_1_flow_rate*self.fluid_spec_heat*inlet_1_T
            if 'inlet_2' in self.inlets_and_outlets['node_'+str(i)]:
                flow_in  = flow_in + inlet_2_flow_rate
                a = a - inlet_2_flow_rate*self.fluid_spec_heat
                b = b + inlet_2_flow_rate*self.fluid_spec_heat*inlet_2_T
            if 'outlet_1' in self.inlets_and_outlets['node_'+str(i)]:
                flow_out = flow_out + inlet_1_flow_rate
            if 'outlet_2' in self.inlets_and_outlets['node_'+str(i)]:
                flow_out = flow_out + inlet_2_flow_rate
            ## Case 1: Top node of the tank
            if i == 1:
                #flow_below = flow coming from the node below (if positive) or going to the node below (if negative)
                flow_below = flow_out - flow_in
                a = a - self.cross_Area*self.top_loss_coeff
                a = a - self.fluid_thermal_conductivity*self.cross_Area/self.node_height
                b = b + self.cross_Area*self.top_loss_coeff*T_out_top
                b = b + self.fluid_thermal_conductivity*self.cross_Area*previous_result['node_'+str(i+1)]['T_ave']/self.node_height
                if flow_below > 0:
                    a = a - flow_below*self.fluid_spec_heat
                    b = b + flow_below*self.fluid_spec_heat*previous_result['node_'+str(i+1)]['T_ave']
                ## flow_above = -flow_below for the next node
                flow_above = -flow_below
            ## Case 2: Bottom node of the tank
            elif i == self.N_nodes:
                a = a - self.cross_Area*self.bottom_loss_coeff
                a = a - self.fluid_thermal_conductivity*self.cross_Area/self.node_height
                b = b + self.cross_Area*self.bottom_loss_coeff*T_out_bottom
                b = b + self.fluid_thermal_conductivity*self.cross_Area*previous_result['node_'+str(i-1)]['T_ave']/self.node_height
                ## since this is the bottom node, a warning is printed if 
                ## the flow coming from or going to a node below is larger than zero
                if abs(flow_above + flow_in - flow_out) > 0.0001:
                    print('Warning: A problem has occurred in the internal flow balance of the tank')
                if flow_above > 0:
                    a = a - flow_above*self.fluid_spec_heat
                    b = b + flow_above*self.fluid_spec_heat*previous_result['node_'+str(i-1)]['T_ave']
            ## Case 3: Intermediate nodes of the tank
            else:
                ## flow_below (flow coming from the node below (if positive) or going to the node below (if negative)) is computed
                ## flow_above has been defined when computing the previous node (node above)
                flow_below = flow_out - flow_in - flow_above
                a = a - 2*self.fluid_thermal_conductivity*self.cross_Area/self.node_height
                b = b + self.fluid_thermal_conductivity*self.cross_Area*previous_result['node_'+str(i+1)]['T_ave']/self.node_height
                b = b + self.fluid_thermal_conductivity*self.cross_Area*previous_result['node_'+str(i-1)]['T_ave']/self.node_height
                if flow_above > 0:
                    a = a - flow_above*self.fluid_spec_heat
                    b = b + flow_above*self.fluid_spec_heat*previous_result['node_'+str(i-1)]['T_ave']
                if flow_below > 0:
                    a = a - flow_below*self.fluid_spec_heat
                    b = b + flow_below*self.fluid_spec_heat*previous_result['node_'+str(i+1)]['T_ave']
                ## flow_above = -flow_below for the next node
                flow_above = -flow_below
            a = a/self.node_thermal_capacity
            b = b/self.node_thermal_capacity
            ## T_final: temperature of the node at the end of the time_step
            ## T_ave: average temperature during the time_step.
            T_final = (self.T['node_'+str(i)]+b/a)*exp(a*time_step)-b/a
            T_ave = (1/(a*time_step))*(self.T['node_'+str(i)]+b/a)*(exp(a*time_step)-1)-b/a
            Sol['node_'+str(i)] = {'T_ave': T_ave, 'T_final': T_final}
        return Sol
    ## Tank.compute outputs: given the inputs of the tank, this function calls the "iteration" function
    ## until convergence is achieved.
    def compute_outputs(self, inlet_1_flow_rate, inlet_1_T, inlet_2_flow_rate, inlet_2_T,
                        T_out_top, T_out_edge, T_out_bottom, time_step, pressure, tolerance = 1e-8):
        ## the tempratures to initialize the iterations are the temperatures of the nodes at the beginning of the time step
        previous_result = {'node_'+str(i): {'T_ave':self.T['node_'+str(i)], 'T_final': self.T['node_'+str(i)]}
                           for i in range(1, self.N_nodes + 1)}
        it = 0
        ## loop to call the iteration function
        while True:
            new_result = self.iteration(inlet_1_flow_rate, inlet_1_T, inlet_2_flow_rate, inlet_2_T, T_out_top,
                                        T_out_edge, T_out_bottom, previous_result, time_step)
            it = it + 1
            ## when convergence is achieved, the loop is ended
            if convergence(previous_result, new_result, tolerance):
                break
            if it == 100:
                print("Maximal number of allowed iterations reached")
                break
            previous_result = new_result
        ## define the outputs
        T_Finals = {}
        for i, node in enumerate(new_result):
            T_Finals[i+1] = new_result[node]['T_final']
        if not temperatures_ordered(T_Finals):
            T_Finals = mix_nodes(T_Finals)
        ## Tank.new_temperatures is a provisional dictionary to save the temperatures that the algorithm determined
        ## this temperatures are set as definitive temperatures of the tank only if the rest of the system has converged as well
        self.new_temperatures = {'node_'+str(i): T_Finals[i] for i in T_Finals}
        for i in range(1, self.N_nodes + 1):
            if 'outlet_1' in self.inlets_and_outlets['node_'+str(i)]:
                outlet_1_T = T_Finals[i]
            if 'outlet_2' in self.inlets_and_outlets['node_'+str(i)]:
                outlet_2_T = T_Finals[i]
        list_average = [ T_Finals[i] for i in T_Finals ]
        mean_temp = np.mean(list_average)
        outputs = {'Outlet 1 Flow Rate': inlet_1_flow_rate, 'Outlet 1 Temperature': outlet_1_T,
                   'Outlet 2 Flow Rate': inlet_2_flow_rate, 'Outlet 2 Temperature': outlet_2_T,
                   'Mean Tank Temperature': mean_temp, 'P': pressure}
        for i in T_Finals:
            outputs['Node '+str(i)+' Temperature'] = T_Finals[i]
        outputs['Iterations'] = it
        return outputs
    ## update_temperature: only when confirming that the inputs and ouputs of the whole system have converged,
    ## the temperature of the tank is updated to the value that it previously determined
    def update_temperature(self):
        self.T = self.new_temperatures
        



    