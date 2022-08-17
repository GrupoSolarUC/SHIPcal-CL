#!/usr/bin/env python
# coding: utf-8

# In[33]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pvlib as pv
import CoolProp.CoolProp as CP
import warnings
from scipy import interpolate
from simulation_functions import corr_exp_solar, Irradiance_2,State


class Heat_exchanger:
    
    def __init__(self,eff,HX_type=None): 
        
        '''
        Heat exchanger's constructor function
        '''
        self.eff=eff
        self.HX_type=HX_type
        
        
    def State_out_calculus(self,m_hot,State_hot_in,m_cold,State_cold_in,P_loss_hot=0,P_loss_cold=0):
        
        '''
        Function to calculate the fluid's state variables at the hot and cold outlet of 
        a heat exchanger given two mass flows.
        
        Parameters:
        ---------
            -m_hot
            -State_hot_in
            -m_cold
            -State_cold_in
            -P_loss_hot
            -P_loss_cold
        
        Outputs:
        ---------
            -State_hot_out
            -State_cold_out
        '''
        
        T_hot_in=State_hot_in['T']    # Temperature at the hot side's inlet
        T_cold_in=State_cold_in['T']  # Temperature at the cold side's inlet
        
        h_hot_in=State_hot_in['H']    # Temperature at the hot side's inlet
        h_cold_in=State_cold_in['H']  # Temperature at the cold side's inlet
        
        P_hot_out=State_hot_in['P']-P_loss_hot    # Pressure at the hot side's outlet
        P_cold_out=State_cold_in['P']-P_loss_cold # Pressure at the cold side's outlet
        
        cp_hot=State_hot_in['C']   # Specific heat capacity at the hot side's inlet  
        cp_cold=State_cold_in['C'] # Specific heat capacity at the cold side's inlet
        
        fluid_hot=State_hot_in['fluid']   # Fluid at the hot side
        fluid_cold=State_cold_in['fluid'] # Fluid at the cold side

        
        it=0 #iteration counter
        
        while True: # Iteration loop for the specific heat capacity at hot and cold sides
        
            C_hot=cp_hot*m_hot    # Product between specific heat capacity and mass flow of the hot side
            C_cold=cp_cold*m_cold # Product between specific heat capacity and mass flow of the cold side

            C_min=min(C_cold,C_hot) # Minimum value of C

            q=self.eff*C_min*(T_hot_in-T_cold_in) # Transfered heat in the heat exchanger
            
            m_list=[m_hot,m_cold]
            C_list=[C_hot,C_cold]
            q_h=self.eff*m_list[C_list.index(min(C_list))]*(h_hot_in-h_cold_in)

            T_hot_out=T_hot_in-q/C_hot     # Temperature at the hot side's outlet
            T_cold_out= T_cold_in+q/C_cold # Temperature at the cold side's outlet
        
            
            State_hot_out=State(fluid_hot,'P',P_hot_out,'T',T_hot_out)     # Variables of state at the hot side's outlet
            State_cold_out=State(fluid_cold,'P',P_cold_out,'T',T_cold_out) # Variables of state at the hot side's outlet
            

            err_hot=(State_hot_out['C']-cp_hot)/State_hot_out['C']*100     # Porcentual error of the cp used at the hot side
            err_cold=(State_cold_out['C']-cp_cold)/State_cold_out['C']*100 # Porcentual error of the cp used at the cold side
            
            if (np.abs(err_hot)<0.1 and np.abs(err_cold)<0.1) or it>=200:  #Criteria for quit the loop

                break

            else:
                it+=1
                cp_hot=(cp_hot+State_hot_out['C'])/2    # New specific heat capacity at the hot side
                cp_cold=(cp_cold+State_cold_out['C'])/2 # New specific heat capacity at the cold side
        
        return State_hot_out,State_cold_out
        return {'State_hot_out':State_hot_out,'State_cold_out':State_cold_out}


