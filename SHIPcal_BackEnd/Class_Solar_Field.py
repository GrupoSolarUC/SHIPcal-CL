# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 17:08:53 2022

@author: Adrian
"""

from scipy.integrate import dblquad
from scipy.interpolate import interp1d
from math import pi
from numpy import arctan, tan, sin, cos, log, exp
import numpy as np
import sys




def diffuse_iam(IAM, beta):
    if IAM == 0:
        return {'sky': 1, 'ground': 1}
    if IAM < 0:
        def IAM_func(theta, phi):
            transverse_angle = arctan(tan(theta)*sin(phi))
            if transverse_angle >= pi/2:
                return 0
            if transverse_angle > pi/3:
                return (abs(IAM) + 1)*(pi/2 - transverse_angle)/(pi/2 - pi/3)
            return 1 - IAM*(1/cos(transverse_angle) - 1)
    if IAM > 0:
        def IAM_func(theta, phi):
            if theta >= pi/2:
                return 0
            IAM_tentative = 1 - IAM*(1/cos(theta) - 1)
            if IAM_tentative <= 0:
                return 0
            return IAM_tentative
    def function_to_integrate(theta, phi, mode):
        if mode == 'numerator':
            return IAM_func(theta, phi)*cos(theta)*sin(theta)
        if mode == 'denominator':
            return cos(theta)*sin(theta)
    if beta == 0 or beta == 90:
        numerator = dblquad(function_to_integrate, 0, pi/2, 0, pi/2, args = ('numerator',))[0]
        denominator = dblquad(function_to_integrate, 0, pi/2, 0, pi/2, args = ('denominator',))[0]
        return {'sky': numerator/denominator, 'ground': (beta == 90)*numerator/denominator}
    if beta < 0 or beta > 90:
        print('La inclinación de los colectores debe estar en el rango [0° , 90°]')
        sys.exit()
    def integration_limit(phi):
        return arctan(tan(pi/2 - beta*pi/180)/cos(phi))
    numerator_sky_1 = dblquad(function_to_integrate, 0, pi/2, 0, pi/2, args = ('numerator',))[0]
    denominator_sky_1 = dblquad(function_to_integrate, 0, pi/2, 0, pi/2, args = ('denominator',))[0]
    numerator_sky_2 = dblquad(function_to_integrate, 0, pi/2, 0, integration_limit, args = ('numerator',))[0]
    denominator_sky_2 = dblquad(function_to_integrate, 0, pi/2, 0, integration_limit, args = ('denominator',))[0]
    numerator_ground = dblquad(function_to_integrate, 0, pi/2, integration_limit, pi/2, args = ('numerator',))[0]
    denominator_ground = dblquad(function_to_integrate, 0, pi/2, integration_limit, pi/2, args = ('denominator',))[0]
    return {'sky': (numerator_sky_1 + numerator_sky_2)/(denominator_sky_1 + denominator_sky_2),
            'ground': numerator_ground/denominator_ground}
        

class Solar_Field:
    '''
    Model for a thermal solar field. The model considers energy losses which are produced due to vapor production and expuslion (this is considered only when pure water is used as working fluid).
    When vapor is produced due to boiling inside the collectors, a flow of mains water is introduced in the inlet of the corresponding row in order to compensate the loss of mass due to vapor expulsion.
    
    Initialization parameters:
        - coll_L: length of each individual collector (m)
        - coll_W: width of each individual collector (m)
        - coll_n0: optical efficiency of the collectors (no units)
        - coll_a1: Linear thermal loss coefficient of the collectors ( kJ / (hr * m^2 * K) )
        - coll_a2: Quadratic thermal loss coefficient of the collectors ( kJ / (hr * m^2 * K^2) )
        - coll_rows: Number of collector rows in the solar field (rows are arranged in parallel; within each row, the collectors are arranged in series) (no units)
        - colls_per_row: Number of collectors in each row (no units)
        - fluid_props: Object of the class "Properties" to define the fluid that is used as working fluid to extract heat from the collectors.
        - tolerance: Maximum discrepancy between successive iterations which is allowed to yield the outputs of the field (iterations are only needed when vapor is being generated inside the collectors) (no units)
        - maxIt: Maximum number of iterations permitted until the outputs are yielded despite convergence has not been achieved (no units)
        
    Attributes (See initialization parameters for description):
        - self.coll_L
        - self.coll_W
        - self.coll_A = self.coll_L * self.coll_W  #area of the collector
        - self.n0
        - self.a1
        - self.a2
        - self.coll_rows
        - self.colls_per_row
        - self.fluid_props
        - self.tolerance
        - self.maxIt
        
    Methods:
        - self.convergence: evaluates whether convergence has been achieved, in the case that iterations are needed due to vapor generation.
        - self.h_out_coll: output computation for a single collector
        - self.h_out_row: output computation for a row of collectors
        - self.compute_outputs: output computation for the whole field
    '''
    def __init__(self, coll_A, coll_n0, coll_a1, coll_a2, IAM_b0, coll_test_flow_per_m2, test_fluid_props,
                 coll_tilt, coll_azimuth, coll_rows, colls_per_row, fluid_props, tolerance = 1e-6, maxIt = 100):
        self.coll_A = coll_A
        self.n0 = coll_n0
        self.a1 = coll_a1
        self.a2 = coll_a2
        self.coll_rows = coll_rows
        self.colls_per_row = colls_per_row
        self.fluid_props = fluid_props
        phi = coll_azimuth*pi/180
        beta = coll_tilt*pi/180
        self.rotation_matrix = np.linalg.inv(np.array([[cos(phi)*cos(beta), -sin(phi), cos(phi)*sin(beta)],
                                                       [sin(phi)*cos(beta),  cos(phi), sin(phi)*sin(beta)],
                                                       [-sin(beta),          0,        cos(beta)         ] ]))
        diff_iam = diffuse_iam(IAM_b0, coll_tilt)
        self.sky_diffuse_iam = diff_iam['sky']
        self.ground_diffuse_iam = diff_iam['ground']
        self.IAM_b0 = IAM_b0
        test_flow = coll_test_flow_per_m2*coll_A
        self.FpUL = - (test_flow*test_fluid_props.fluid_spec_heat/coll_A)*log(1 - coll_a1*coll_A/(test_flow*test_fluid_props.fluid_spec_heat))
        self.test_conditions = ( test_flow*test_fluid_props.fluid_spec_heat/(coll_A*self.FpUL) )*( 1 - exp(- coll_A*self.FpUL/(test_flow*test_fluid_props.fluid_spec_heat) ))
        self.coll_azimuth = phi
        self.coll_tilt = beta
        self.tolerance = tolerance
        self.maxIt = maxIt
    def convergence(self, m_mains, m_vap):
        """
        Function that evaluates the convergence of the iterations (iterations are only necessary when vapor is being generated).
        In each iteration, a flow of mains water is introduced into the collector row to compensate the mass loss due to vapor expulsion.
        In order to achieve convergence, the mass flow of water being introduced and the mass flow of vapor being produced must be equal (or have a discrepancy smaller than a "tolerance" parameter)

        Parameters:
            - m_mains: mass flow of mains water that was introduced into the collector at the beginning of the iteration (kg/hr).
            - m_vap: mass flow of vapor produced due to the heat of the collectors (kg/hr).

        Returns:
            a boolean value which is True if convergence has been achieved, and False otherwise.
        """
        if m_mains == m_vap:
            return True
        m_max = max([ m_mains, m_vap ])
        return abs(m_mains - m_vap)/m_max < self.tolerance
    def IAM_func(self, angle_of_incidence, longitudinal_angle, transverse_angle):
        if self.IAM_b0 == 0:
            return 1
        if self.IAM_b0 < 0:
            if transverse_angle >= pi/2:
                return 0
            if transverse_angle > pi/3:
                return (abs(self.IAM_b0) + 1)*(pi/2 - transverse_angle)/(pi/2 - pi/3)
            return 1 - self.IAM_b0*(1/cos(transverse_angle) - 1)
        if self.IAM_b0 > 0:
            if angle_of_incidence >= pi/2:
                return 0
            IAM_tentative = 1 - self.IAM_b0*(1/cos(angle_of_incidence) - 1)
            if IAM_tentative <= 0:
                return 0
            return IAM_tentative
    def incidence_angles(self, solar_zenith, solar_azimuth):
        if solar_zenith == 0:
            solar_N = 0
            solar_E = 0
            solar_Z = 1
        elif solar_zenith == 180:
            solar_N = 0
            solar_E = 0
            solar_Z = -1
        else:
            solar_N = cos(solar_azimuth*pi/180)*sin(solar_zenith*pi/180)
            solar_E = sin(solar_azimuth*pi/180)*sin(solar_zenith*pi/180)
            solar_Z = cos(solar_zenith*pi/180)
        coll_coordinates = np.dot(self.rotation_matrix, np.array([[solar_N],
                                                                  [solar_E],
                                                                  [solar_Z]]))
        coll_N = coll_coordinates[0][0]
        coll_E = coll_coordinates[1][0]
        coll_Z = coll_coordinates[2][0]
        if coll_Z <= 0:
            return None, None, None
        if coll_E == 0 and coll_N == 0:
            trans = 0
            longi = 0
            aoi = 0
        elif coll_E == 0:
            trans = 0
            longi = arctan(abs(coll_N)/coll_Z)
            aoi = longi
        elif coll_N == 0:
            trans = arctan(abs(coll_E)/coll_Z)
            longi = 0
            aoi = trans
        else:
            trans = arctan(abs(coll_E)/coll_Z)
            longi = arctan(abs(coll_N)/coll_Z)
            tg2_aoi = (tan(trans))**2 + (tan(longi))**2
            aoi = arctan(tg2_aoi**(1/2))
        if longi > pi/2 - self.coll_tilt and coll_N > 0:
            longi = pi/2 - self.coll_tilt
            tg2_aoi = (tan(trans))**2 + (tan(longi))**2
            aoi = arctan(tg2_aoi**(1/2))
        return aoi, longi, trans
    def h_out_coll(self, m_in, h_in, rad, IAM_eff, T_amb, flow_correction_factor):
        '''
        Function that computes the specific enthalpy of the flow leaving an individual collector.
        It considers the flow temperature and the ambient temperature in order to estimate the efficiency of the collector.

        Parameters:
            - m_in: mass flow entering the collector (kg/hr)
            - h_in: specific enthalpy of the flow entering the collector (kJ/kg)
            - rad: solar radiation reaching the collector ( kJ / ( hr * m^2 ) )
            - T_amb: ambient temperature (°C)

        Returns:
            - specific enthalpy of the flow leaving the collector.
        '''
        T_in = self.fluid_props.h_to_T(h_in)
        if T_amb > T_in:
            eff = self.n0*IAM_eff
        else:
            eff = self.n0*IAM_eff - self.a1*(T_in - T_amb)/rad - self.a2*((T_in - T_amb)**2)/rad
        if eff < 0:
            eff = 0
        eff = flow_correction_factor*eff
        power = eff*rad*self.coll_A
        return h_in + power/m_in
    def h_out_row(self, m_in_row, h_in_row, h_mains, rad, IAM_eff, T_amb, flow_correction_factor):
        '''
        Function that computes the specific enthalpy of the flow leaving a row of collectors.

        Parameters:
            - m_in_row: mass flow entering the row (kg/hr)
            - h_in_row: specific enthalpy of the flow entering the row (kJ/kg)
            - h_mains: specific enthalpy of the mains flow that will be introduced into the field in the case that vapor is generated inside the collectors (kJ/kg).
            - rad: solar radiation reaching the collector row ( kJ / ( hr * m^2 ) )
            - T_amb: ambient temperature (°C)

        Returns a tuple with the following elements:
            - specific enthalpy of the flow leaving the row (kJ/kg)
            - power being "wasted" (only when vapor is generated) to heat the mains flow that is introduced to compensate the vapor loss. This flow is heated from its initial enthalpy to the specific enthalpy of vapor at saturation (kJ/hr).
            - vapor mass flow being released (or mains water flow being introduced to compensate the vapor loss) (kg/hr)
            - number of iterations executed to achieve convergence (no units)
            - a boolean value which is True if the maximal number of allowed iterations was reached without achieving convergence, and False otherwise.
        '''
        m_mains = 0
        it = 0
        while True:
            m_in = m_in_row + m_mains
            h_in = (m_in_row*h_in_row + m_mains*h_mains)/m_in
            for i in range(self.colls_per_row):
                h_in = self.h_out_coll(m_in, h_in, rad, IAM_eff, T_amb, flow_correction_factor)
            h_out = h_in
            m_vap = self.fluid_props.quality(h_out)*m_in
            it = it+1
            if self.convergence(m_mains, m_vap):
                maxIt = False
                break
            if it == self.maxIt:
                maxIt = True
                break
            Q_useful = m_in_row*(self.fluid_props.h_sat_liq - h_in_row)
            Q_actual = m_in*(h_out - (m_in_row*h_in_row + m_mains*h_mains)/m_in)
            Q_waste = Q_actual - Q_useful
            m_mains = Q_waste/(self.fluid_props.h_sat_vap - h_mains)
        if m_vap > 0:
            return self.fluid_props.h_sat_liq, m_vap*(self.fluid_props.h_sat_vap - h_mains), m_vap, it, maxIt
        else:
            return h_out, 0, 0, it, maxIt
    def compute_outputs(self, m_in_field, h_in_field, h_mains, rad, IAM_eff, T_amb):
        '''
        Function that computes the specific enthalpy of the flow that is leaving the solar field.

        Parameters:
            - m_in_field: mass flow entering the solar field (kg/hr).
            - h_in_field: specific enthalpy of the flow entering the solar field (kJ/kg).
            - h_mains: specific enthalpy of the mains flow that will be introduced into the field in the case that vapor is generated inside the collectors (kJ/kg).
            - rad: solar radiation reaching the solar field ( kJ / ( hr * m^2 ) )
            - T_amb: ambient temperature (°C)

        Returns a dictionary with the following keys:
            - 'h_out': specific enthalpy of the flow leaving the tank (kJ/kg)
            - 'Q_useful': useful energy delivered to the flow (i.e. all energy which is not used to generate vapor) (kJ/hr)
            - 'Q_waste': wasted energy which is used to generate vapor inside the collectors (kJ/hr)
            - 'm_vap': mass flow of vapor being released from the field. (kg/hr)
            - 'iterations': number of iterations that were necessary to achieve convergence between the introduced mains flow and the vapor flow being released (no units)
            - 'maxIt': boolean value which is True if the maximum number of allowed iterations was reached without achieving convergence.
        '''
        if rad <= 0 or m_in_field <= 0:
            return {'h_out': h_in_field,
                    'm_out': m_in_field,
                    'Q_useful': 0,
                    'Q_waste': 0,
                    'm_vap': 0,
                    'iterations': 0,
                    'IAM': IAM_eff,
                    'maxIt': False }
        m_in_rows = m_in_field/self.coll_rows
        use_conditions = ( m_in_rows*self.fluid_props.fluid_spec_heat/(self.coll_A*self.FpUL) )*( 1 - exp( - self.coll_A*self.FpUL/( m_in_rows*self.fluid_props.fluid_spec_heat ) ) )
        flow_correction_factor = use_conditions/self.test_conditions
        h_out, Q_waste_rows, m_waste_rows, it, maxIt = self.h_out_row(m_in_rows, h_in_field, h_mains, rad, IAM_eff, T_amb, flow_correction_factor)
        Q_waste = Q_waste_rows*self.coll_rows
        Q_useful = m_in_field*(h_out - h_in_field)
        m_waste = m_waste_rows*self.coll_rows
        return {'h_out': h_out,
                'm_out': m_in_field,
                'Q_useful': Q_useful,
                'Q_waste': Q_waste,
                'm_vap': m_waste,
                'iterations': it,
                'IAM': IAM_eff,
                'maxIt': maxIt}
