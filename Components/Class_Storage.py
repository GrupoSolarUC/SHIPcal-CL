#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pvlib as pv
import CoolProp.CoolProp as CP
import warnings
from scipy import interpolate
from simulation_functions import corr_exp_solar, Irradiance_2,State
from bokeh.plotting.figure import Figure
from bokeh.plotting import figure, output_notebook, show
from bokeh.palettes import Dark2_5 as palette
# itertools handles the cycling
import itertools  


class Storage:
    
    def __init__(self,Storage_param):
        
        '''
        Storage's constructor function
        '''
        self.n_nodes=Storage_param['n_nodes'] #Number of nodes in the storage
        self.T_max=Storage_param['T_max']     #Max temperature of the start temperature profile, °C
        self.T_min=Storage_param['T_min']     #Min temperature of the start temperature profile, °C
        
        self.nodes_Temperature=np.linspace(self.T_max,self.T_min,self.n_nodes) #Start temperature profile
        
        self.h_nodes=Storage_param['H']/self.n_nodes #Lenght of each node, m
        self.P=Storage_param['P']*100000 # Tank´s pressure, Pa
        self.D=D #Diamameter of the tank, m
        self.H=H #Height of the tank, m
        self.A_circ=np.pi*(self.D/2)**2 #Tank's transversal area, m^2
        self.A_cil_ext=np.pi*self.D*self.h_nodes #Cilinder area of each node, m^2
        
        self.k_storage=Storage_param['k_storage'] #Storage´s thermal conductivity, W/mK
        self.e_storage=Storage_param['e_storage'] #Storage´s thickness
        self.k_insulation=Storage_param['k_insulation'] #Insulation´s thermal conductivity, W/mK
        self.e_insulation=Storage_param['e_insulation'] #Insulation´s thickness
        
        self.nodes_height=np.linspace(0,self.H,self.n_nodes+1) #Height of each node, m
        
        self.h_out=Storage_param['h_out'] #Storage´s outlet height, m
        
        
    def M_vector(self,T_in):
        '''
        Funtion that calculates in which node there will be located the inlets and outlets of the tank.
        To the inlets, the node will be determinated by the temperature, defining that the mass flow will enter
        in the node whose temperature is nearest to the temperature of the mass flow. 
        For the outlets, the node will be determined by the height of the outlets define at the moment of definig
        the storage. 
        
        
        
        Parameters:
        -----------
            T_in: Temperature of the inlets flux, °C
            
        Outputs:
        -----------
            self.M_in : Matrix with rows equal to the number of nodes and columns equal to the number of inlets.
                        The location of inlets are marked with a number 1.
            self.M_out: Matrix with rows equal to the number of nodes and columns equal to the number of outlets.
                        The location of inlets are marked with a number -1.
        '''
        
        m_in=np.zeros([self.n_nodes,len(T_in)])
        m_out=np.zeros([self.n_nodes,len(self.h_out)])
        
        for i in range(len(T_in)):
            for j in range(self.n_nodes-1):
                if self.nodes_Temperature[j]>=T_in[i]>self.nodes_Temperature[j+1]:
                    m_in[j][i]=1
                elif self.nodes_Temperature[0]<T_in[i]:
                    m_in[0][i]=1
                elif self.nodes_Temperature[-1]>=T_in[i]:
                    m_in[-1][i]=1
                    
        
        for i in range(len(self.h_out)):
            for j in range(self.n_nodes):
                if self.nodes_height[j]<=self.h_out[i]<self.nodes_height[j+1]:
                    m_out[j][i]=-1
        
        
        self.M_in,self.M_out=m_in,np.flip(m_out)
        
    def R_k_cil_external(self,k,D_ext,D_int):
        
        '''
        Method that calculates the conductive thermal resistance of a cylindrical shell
        
        Parameters:
        ----------- 
            -k    : Thermal conductivity of the cylindrical shell material, W/(m*K)
            -D_ext: External diameter of the cylindrical shell, m
            -D_int: Internal diameter of the cylindrical shell, m
            
        Outputs:
        -----------
            -R_k : Conductive thermal resistance ot the cylindrical shell, K/W
            
        '''
        L=self.h_nodes # Height of the cylindrical shell, m
        
        R_k=np.log(D_ext/D_int)/(2*np.pi*L*k) #Conductive thermal resistance ot the cylindrical shell, K/W
        return R_k
    
    def UA_top_bottom(self,T_inf,T_f,g=9.8):
        
        '''
        Function that calculates the UA between the fluid inside the tank and the environment at the bottom and the upper tank's faces.
        Assuming that the internal face of the storage has the same temperature that the fluid. 
        
        Parameters:
        -----------
            T_inf: Environmental temperature, °C
            T_f : Fluid's temperature, °C
            
        Outputs:
        -----------
            UA: UA between the fluid inside the tank and the environment, W/K
        '''
        
        L=self.h_nodes
        D=self.D
        
        T_s=(T_f+T_inf)/2
        
        it=0
        while True:    
            
            T_p=(T_inf+T_s)/2
            nu_air_ext=CP.PropsSI('V','P',self.P,'T',T_p+273,'air')/CP.PropsSI('D','P',self.P,'T',T_p+273,'air') 
            B_air_ext=1/(T_p+273) 
            Gr=g*B_air_ext*np.abs(T_s-T_inf)*(L**3)/(nu_air_ext**2)
            Pr_air_ext=CP.PropsSI('PRANDTL','P',self.P,'T',T_p+273,'air')
            Ra=Gr*Pr_air_ext

            if 0<Ra<10**9:#10**4<Ra<10**9:
                h=1.32*(np.abs(T_s-T_inf)/D)**0.25
            elif Ra>=10**9:
                h=1.43*np.abs(T_s-T_inf)**(1/3)
            else:
                h=0.001
                
            if self.k_storage!=0: 
                R_k_storage=self.e_storage/(self.k_storage*self.A_circ)
            else:
                R_k_storage=0
            
            if self.k_insulation!=0: 
                R_k_insulation=self.e_insulation/(self.k_insulation*self.A_circ)    
            else:
                R_k_insulation=0
                
            UA=(1/(h*self.A_cil_ext)+R_k_storage+R_k_insulation)**(-1)
                  
            Q_f_inf=(T_f-T_inf)*UA
            Q_s_inf=(T_s-T_inf)*self.A_circ*h
            
            err=np.abs(Q_f_inf-Q_s_inf)/Q_s_inf*100
            
            if err<=0.1 or it>=200:
                break
            
            else:
                T_s=T_inf+Q_f_inf/(self.A_circ*h)
                it+=1
                
        return UA
        
        
        
    def UA_external_surface(self,T_inf,T_f,g=9.8):
        '''
        Function that calculates the UA between the fluid inside the tank and the environment at the external cilinder's face.
        Assuming that the internal face of the storage has the same temperature that the fluid. 
        
        Parameters:
        -----------
            T_inf: Environmental temperature, °C
            T_f : Fluid's temperature, °C
            
        Outputs:
        -----------
            UA_ext: UA between the fluid inside the tank and the environment, W/K
        '''
        L=self.h_nodes
        D=self.D
        
        T_s=(T_f+T_inf)/2
        
        it=0

        while True:    
            
            T_p=(T_inf+T_s)/2
            nu_air_ext=CP.PropsSI('V','P',self.P,'T',T_p+273,'air')/CP.PropsSI('D','P',self.P,'T',T_p+273,'air')
            B_air_ext=1/(T_p+273)   
            Gr=g*B_air_ext*np.abs(T_s-T_inf)*(L**3)/(nu_air_ext**2)
            Pr_air_ext=CP.PropsSI('PRANDTL','P',self.P,'T',T_p+273,'air')
            Ra=Gr*Pr_air_ext

            if 10**4<Ra<10**9:#10**4<Ra<10**9:
                h=1.42*(np.abs(T_s-T_inf)/L)**0.25
            elif Ra>=10**9:
                h=0.95*np.abs(T_s-T_inf)**(1/3)
            else:
                h=0.001
               
            if self.k_storage!=0: 
                R_k_storage=self.R_k_cil_external(self.k_storage,D_ext=D,D_int=D-self.e_storage)    
            else:
                R_k_storage=0
            
            if self.k_insulation!=0: 
                R_k_insulation=self.R_k_cil_external(self.k_insulation,D_ext=D+self.e_insulation,D_int=D)    
            else:
                R_k_insulation=0
                    
            UA_ext=(1/(h*self.A_cil_ext)+R_k_storage+R_k_insulation)**(-1)
            
            Q_f_inf=(T_f-T_inf)*UA_ext
            Q_s_inf=(T_s-T_inf)*self.A_cil_ext*h
            
            err=np.abs((Q_f_inf-Q_s_inf)/Q_s_inf*100)
            if err<=0.1 or it>=200:
                break
            
            else:
                T_s=T_inf+Q_f_inf/(self.A_cil_ext*h)
                it+=1
                
        return UA_ext
    
    def Conduction(self,T_hot,T_cold):
        '''
        Function that calculates the heat transferred from one node to another by conduction.
        
        Parameters:
        ----------
            T_hot: Temperature of the hot node, °C
            T_cold: Temperature of the cold node, °C
        
        Outputs:
        --------
            Q_cond: heat transferred from one node to another by conduction, W
        '''
        
        k_nodes=CP.PropsSI('L','P',self.P,'T',np.linspace(273+np.max(self.nodes_Temperature),273+np.min(self.nodes_Temperature),n_nodes),'water')
        k_nodes=np.mean(k_nodes)
        
        Q_cond=k_nodes*self.A_circ*(T_hot-T_cold)/self.h_nodes
        
        return Q_cond
    
    def Mass_balance(self,m_in,m_out):
        '''
        Function that calculates one of the inlets mass flow to enforce the mass conservation inside the tank. This function assumes that 
        there are the exact same number of inlets and outlets. 
        
        Parameters:
        -----------
        m_in: Array or list that contains the values of the inlets mass flow.
        m_out: Array or list that contains the values of the outlets mass flow.
        '''
        
        M_in=np.sum(m_out)-np.sum(m_in)
        m_in=np.insert(m_in,len(m_in),M_in)
        return(m_in)
    
    def Energy_balance(self,dt,T_inf):
        '''
        Function that calculates and updates the temperature of each node by doing an energy balance.
        
        Parameters:
        -----------
        dt   : Time step of the simulation, s.
        T_inf: Ambient temperature, °C
        
        Outputs:
        ---------
        nodes.Temperature: Numpy array of shape (n_nodes,) that contains the resulting temperature of each node without 
        considering the mass flow entering and leaving the tank.
        
        '''        
        
        cp_nodes=CP.PropsSI('C','P',self.P,'T',                               # Array that contains 
                            np.linspace(273+np.max(self.nodes_Temperature),   # the value of Cp on
                            273+np.min(self.nodes_Temperature),n_nodes),      # each node
                            'water')                                          #
        cp_nodes=np.mean(cp_nodes) # Mean value of Cp, J/(kg*K)
        
        rho_nodes=CP.PropsSI('D','P',self.P,'T',                              # Array that contains
                             np.linspace(273+np.max(self.nodes_Temperature),  # the value of density 
                             273+np.min(self.nodes_Temperature),n_nodes),     # on each node
                             'water')                                         #
        rho_nodes=np.mean(rho_nodes) # Mean value of density, kg/m^{3}
        
        V=self.A_circ*self.h_nodes  #Volume of each node, m^{3}
         
        new_nodes_Temperature=np.zeros(n_nodes) #Empty array to fulfill with the new temperature of the nodes 
        
        for i in range(self.n_nodes): # Loop to go through all the nodes
            
            if i>=1 and i<n_nodes-1:  #Energy balance for the middle nodes
                
                Q_cond=(self.Conduction(self.nodes_Temperature[i-1],self.nodes_Temperature[i])   # Heat transferred to the current node from
                                                                                                 # the node located above it
                        
                        -self.Conduction(self.nodes_Temperature[i],self.nodes_Temperature[i+1])) # Heat transferred from the current node 
                                                                                                 # to the node located below it
                    
                Q_loss_external_surface=(-self.UA_external_surface(T_inf,self.nodes_Temperature[i]) # Heat transferred from the storage´s inside
                                         *(self.nodes_Temperature[i]-T_inf))                        # to the environment

                
                new_nodes_Temperature[i]=self.nodes_Temperature[i]+(Q_cond+Q_loss_external_surface)*dt/(V*rho_nodes*cp_nodes) # New temperature of the node i
                 
            elif i==0: #Energy balance for the node at the top of the storage     
                
                Q_cond=-self.Conduction(self.nodes_Temperature[i],self.nodes_Temperature[i+1]) # Heat transferred from the current node
                                                                                               # to the node located below it
                    
                Q_loss_external_surface=(-self.UA_external_surface(T_inf,self.nodes_Temperature[i]) # Heat transferred from the storage´s inside
                                         *(self.nodes_Temperature[i]-T_inf))                        # to the environment by the cylindrical surface
                
                Q_loss_top=(-self.UA_top_bottom(T_inf,self.nodes_Temperature[i]) # Heat transferred from the storage´s inside
                            *(self.nodes_Temperature[i]-T_inf))                  # to the environment by the top surface
                
                
                new_nodes_Temperature[i]=self.nodes_Temperature[i]+(Q_cond+Q_loss_external_surface+Q_loss_top)*dt/(V*rho_nodes*cp_nodes) # New temperature of the node i
                
            else: #Energy balance for the node at the top of the storage     
                
                Q_cond=self.Conduction(self.nodes_Temperature[i-1],self.nodes_Temperature[i]) # Heat transferred to the current node from
                                                                                              # the node located above it
                
                Q_loss=(-self.UA_external_surface(T_inf,self.nodes_Temperature[i]) # Heat transferred from the storage´s inside
                        *(self.nodes_Temperature[i]-T_inf))                        # to the environment by the cylindrical surface
                                                                                                                    
                Q_loss_bottom=(self.UA_top_bottom(T_inf,self.nodes_Temperature[i]) # Heat transferred between the node at the bottom 
                               *(self.nodes_Temperature[i]-T_inf))                 # of the storage and the enviroment
                
                new_nodes_Temperature[i]=self.nodes_Temperature[i]+(Q_cond+Q_loss_external_surface+Q_loss_bottom)*dt/(V*rho_nodes*cp_nodes) # New temperature of the node i
                
        self.nodes_Temperature=new_nodes_Temperature # Update of the nodes temperature
        return self.nodes_Temperature
            
    def Energy_balance_mass_flow(self,dt,T_inf):
        '''
        Function that calculates and updates the temperature of each node by doing an energy balance, considering the mass flow entering and 
        leaving the tank.
        
        Parameters:
        -----------
        dt   : Time step of the simulation, s.
        T_inf: Ambient temperature, °C
        
        Outputs:
        ---------
        nodes.Temperature: Numpy array of shape (n_nodes,) that contains the resulting temperature of each node
        
        '''        
        
        cp_nodes=CP.PropsSI('C','P',self.P,'T',                               # Array that contains 
                            np.linspace(273+np.max(self.nodes_Temperature),   # the value of Cp on
                            273+np.min(self.nodes_Temperature),n_nodes),      # each node
                            'water')                                          #
        cp_nodes=np.mean(cp_nodes) # Mean value of Cp, J/(kg*K)
        
        rho_nodes=CP.PropsSI('D','P',self.P,'T',                              # Array that contains
                             np.linspace(273+np.max(self.nodes_Temperature),  # the value of density 
                             273+np.min(self.nodes_Temperature),n_nodes),     # on each node
                             'water')                                         #
        rho_nodes=np.mean(rho_nodes) # Mean value of density, kg/m^{3}
        
        V=self.A_circ*self.h_nodes  #Volume of each node, m^{3}
         
        new_nodes_Temperature=np.zeros(n_nodes) #Empty array to fulfill with the new temperature of the nodes 
        
        for i in range(self.n_nodes): # Loop to go through all the nodes
            
            if i>=1 and i<n_nodes-1:  #Energy balance for the middle nodes
                
                Q_cond=(self.Conduction(self.nodes_Temperature[i-1],self.nodes_Temperature[i])   # Heat transferred to the current node from
                                                                                                 # the node located above it
                        
                        -self.Conduction(self.nodes_Temperature[i],self.nodes_Temperature[i+1])) # Heat transferred from the current node 
                                                                                                 # to the node located below it
                    
                Q_loss_external_surface=(-self.UA_external_surface(T_inf,self.nodes_Temperature[i]) # Heat transferred from the storage´s inside
                                         *(self.nodes_Temperature[i]-T_inf))                        # to the environment

                
                new_nodes_Temperature[i]=self.nodes_Temperature[i]+(Q_cond+Q_loss_external_surface)*dt/(V*rho_nodes*cp_nodes) # New temperature of the node i
                 
            elif i==0: #Energy balance for the node at the top of the storage     
                
                Q_cond=-self.Conduction(self.nodes_Temperature[i],self.nodes_Temperature[i+1]) # Heat transferred from the current node
                                                                                               # to the node located below it
                    
                Q_loss_external_surface=(-self.UA_external_surface(T_inf,self.nodes_Temperature[i]) # Heat transferred from the storage´s inside
                                         *(self.nodes_Temperature[i]-T_inf))                        # to the environment by the cylindrical surface
                
                Q_loss_top=(-self.UA_top_bottom(T_inf,self.nodes_Temperature[i]) # Heat transferred from the storage´s inside
                            *(self.nodes_Temperature[i]-T_inf))                  # to the environment by the top surface
                
                
                new_nodes_Temperature[i]=self.nodes_Temperature[i]+(Q_cond+Q_loss_external_surface+Q_loss_top)*dt/(V*rho_nodes*cp_nodes) # New temperature of the node i
                
            else: #Energy balance for the node at the top of the storage     
                
                Q_cond=self.Conduction(self.nodes_Temperature[i-1],self.nodes_Temperature[i]) # Heat transferred to the current node from
                                                                                              # the node located above it
                
                Q_loss=(-self.UA_external_surface(T_inf,self.nodes_Temperature[i]) # Heat transferred from the storage´s inside
                        *(self.nodes_Temperature[i]-T_inf))                        # to the environment by the cylindrical surface
                                                                                                                    
                Q_loss_bottom=(self.UA_top_bottom(T_inf,self.nodes_Temperature[i]) # Heat transferred between the node at the bottom 
                               *(self.nodes_Temperature[i]-T_inf))                 # of the storage and the enviroment
                
                new_nodes_Temperature[i]=self.nodes_Temperature[i]+(Q_cond+Q_loss_external_surface+Q_loss_bottom)*dt/(V*rho_nodes*cp_nodes) # New temperature of the node i
                
        self.nodes_Temperature=new_nodes_Temperature # Update of the nodes temperature
        return self.nodes_Temperature