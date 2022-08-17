#!/usr/bin/env python
# coding: utf-8

# In[14]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pvlib as pv
import CoolProp.CoolProp as CP
import warnings
from scipy import interpolate
from simulation_functions import corr_exp_solar, Irradiance_2,State


class Collector():
    'Class '
    
    def __init__(self,coll_parameters):
        
        '''Collector´s constructor function'''
        
        self.L=coll_parameters['L']                      #Collector´s length, m
        self.W=coll_parameters['W']                      #Collector´s width, m 
        self.n0=coll_parameters['n0']                    #Collector's optical efficiency
        self.a1=coll_parameters['a1']                    #Collector's linear thermal loss coefficient, W/(m^2*K)
        self.a2=coll_parameters['a2']                    #Collector's quadratic thermal loss coefficient,  W/(m^2*K^2)
        self.price=coll_parameters['price']              #Collector´s price, $CLP
        self.A=coll_parameters['L']*coll_parameters['W'] #Collector´s aperture area, m^2
        
    def n_CP(self,T,T_amb,G=1000):
        '''
        Function that calculates the efficiency of a solar thermal collector considering the losses to the ambient.

        Inputs: 
        -----------
            -T : Fluid´s inlet temperature to the collector, °C
            -Tm: Ambien temperature, °C
            -G : Solar radiation, W/m2

        Output:
        ----------
            -n: Collector's efficiency
        '''

        if G==0:
            return 0
        else:
            n=self.n0-(T-T_amb)*self.a1/G-self.a2*((T-T_amb)**2)/G
            return n
    
    def Salto_Temp(self,m_col,G,cp,T,T_amb):
        '''
        Funtion that calculates the temperature difference between the outlet and the inlet of the collector.

        -Inputs:
        ------------
            -m   :Mass flow of fluid entering the collector, kg/s
            -eta :Collector's efficiency
            -G   :Irradiance at tilted collector surface area, W/m2
            -A   :Collector´s aperture area, m2
            -cp  :Fluid´s specific heat at constat pressure, kJ/(kg*K)

        -Outputs:
        ------------
            -Salto_T: Temperature difference between the outlet and the inlet of the collector , °C
        '''
        eta=self.n_CP(T,T_amb)
        Salto_T=G*eta*self.A/m_col/cp/1000
        return Salto_T


def Solar_field_construction(N_rows,N_col,Coll_parameters):
    
    '''
    Funtion that define a solar field based on the number of rows connected
    in parallel and the number of collector´s conected in series in each row
    
    Inputs:
    --------
        -N_rows         : Number of rows conected in parallel in the solar field
        -N_col          : Number of collectors that are conected in series in each row
        -Coll_parameters: Dictionary that contains the parameters of one single collector:{'L':Lenght [m],'W':width [m],'n0':Optical efficiency,
                                                                                            'a1':Linear thermal loss coefficient [W/(m^2*K)],
                                                                                            'a2':quadratic thermal loss coefficient [W/(m^2*K^2)],
                                                                                            'price':price}     
                                                                                            
    Outputs:
    ---------
        solar_field: Dictionary that contains all the collectors and the number of rows and columns
    '''
    
    solar_field={'N_rows':N_rows,'N_col':N_col} #Creates a dictionary to fill with collectors
    
    for row in range(1,N_rows+1):    # Loop to go trowh the number of rows 
        for col in range(1,N_col+1): # Loop to go trowh the number of columns
            
            solar_field[f'Coll_{row}_{col}']=Collector(Coll_parameters) #Creates each collector wit the name: 
                                                                        #"Coll_{row number}_{column number}"
            
    return solar_field

def Solar_field_operation(solar_field,m_field,state_in,hour,Climate_Data):
    '''
    Function that calculates the output state of a fluid that enters a solar field. 
    
    Inputs:
    ------
        -solar_field : Dictionary that contains all the collectors of the solar field and the numer of rows and columns
        -m_field     : Mass flow of fluid entering the solar field, kg/s
        -state_in    : Dictionary that contains fluid's state variables entering the collector field. This dictanary
                       can be created using the function State().
        -hour        : Hour of the year that the simulation is being perform.
        -Climate_Data: Pandas dataframe that contains the climate file parameters: 'Year', 'Month', 'Day', 'Hour', 
                       'Minute', 'GHI', 'DNI', 'DHI', 'Tdry','Tdew', 'RH', 'Pres', 'Wspd', 'Wdir', 'Snow Depth'.
                      
    Outputs:
    ------
        -state_out: Dictionary that contains the fluid´s state variables at the outlet of the solar field
    '''
    
    N_rows=solar_field['N_rows']
    N_col=solar_field['N_col']
    m_rows=m_field/N_rows #Uniformly distributed the mass flow of the field in each row

    h_out_rows=np.zeros(N_rows) #Empty array for the outlet  specific enthalpy of each row 
    P_field=state_in['P'] #Saves the pressure of the field 
    
    for row in range(1,N_rows+1):    # Loop to go trowh the number of rows 
        for col in range(1,N_col+1): # Loop to go trowh the number of columns
            
            if row!=1:  #Makes the distition if the row is not the first one (the most north located). This is to consider
                        #the radiation losses caused by the self shading of the field. 
                
                T_col_in=state_in['T']                                #Fluid's temperature at the collector's inlet
                cp=state_in['C']                                      #Fluid's specific heat at the collector's inlet
                G=Climate_Data['DNI'][hour]+Climate_Data['DHI'][hour] #Value of the solar radiation at the collector's tilted surface,
                                                                      #take in consideration the self shading caused by the other rows
                T_amb=Climate_Data['Tdry'][hour]                      #Dry bulb temperature 
                
                Delta_T=solar_field[f'Coll_{row}_{col}'].Salto_Temp(m_col=m_rows,G=G, #Calculus of the temperature
                                                                    cp=cp,T=T_col_in, #difference between the inlet
                                                                    T_amb=T_amb)      #and the outlet of the solar field
                
                T_col_out=T_col_in+Delta_T            #Fluid's temperature a the collector's outlet
                state_in=State(state_in['fluid'],'T', #Fludi's state at the collector's outlet
                               T_col_out,'P',P_field) #
                
            else:
                T_col_in=state_in['T']                                #Fluid's temperature at the collector's inlet
                cp=state_in['C']                                      #Fluid's specific heat at the collector's inlet
                G=Climate_Data['DNI'][hour]+Climate_Data['DHI'][hour] #Value of the solar radiation at the collector's tilted surface
                T_amb=Climate_Data['Tdry'][hour]                      #Dry bulb temperature
                
                Delta_T=solar_field[f'Coll_{row}_{col}'].Salto_Temp(m_col=m_rows,G=G,  #Calculus of the temperature
                                                                    cp=cp,T=T_col_in,  #difference between the inlet
                                                                    T_amb=T_amb)       #and the outlet of the solar field
                
                T_col_out=T_col_in+Delta_T            #Fluid's temperature a the collector's outlet
                state_in=State(state_in['fluid'],'T', #Fludi's state at the collector's outlet
                               T_col_out,'P',P_field) #
        
        h_out_rows[row-1]=state_in['H'] #Array that contains the fluid's speficic enthalpy at the end of each row 
    
    h_out_field=np.sum(m_rows*h_out_rows)/(m_field)
    state_out=State(state_in['fluid'],'H',h_out_field,'P',P_field)
    return state_out
    