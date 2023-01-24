# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 09:22:22 2023

@author: adria
"""

from Run_Simulation import simulate_system
from simulation_functions import Heat_Map

system_params = {'sector': 'Ninguno',
                  'sim_name': 'Simulacion',
                  'application': 'Ninguna',
                  'processes': 'Ninguno',
                  'location': 'Santiago',
                  'fuel_name': 'Diesel',
                  'yearly_demand_unit': 'm3',
                  'demand_monday': True,
                  'demand_tuesday': True,
                  'demand_wednesday': True,
                  'demand_thursday': True,
                  'demand_friday': True,
                  'demand_saturday': False,
                  'demand_sunday': False,
                  'operation_start': '8:00',
                  'operation_end': '18:00',
                  'demand_january': 1.5,
                  'demand_february': 1.5,
                  'demand_march': 1.5,
                  'demand_april': 1.5,
                  'demand_may': 1.8,
                  'demand_june': 1.8,
                  'demand_july': 1.8,
                  'demand_august': 1.8,
                  'demand_september': 1.5,
                  'demand_october': 1.5,
                  'demand_november': 1.5,
                  'demand_december': 1.5,
                  'boiler_nominal_power': 500,
                  'boiler_nominal_power_units': 'kW',
                  'boiler_pressure': 5,
                  'boiler_pressure_units': 'bar',
                  'boiler_type': 'Condensación líquido',
                  'boiler_efficiency': 0.8,
                  'closed_system': False,
                  'monthly_return_temperature': None,
                  'return_inlet_temperature': None,
                  'temperature_january': None,
                  'temperature_february': None,
                  'temperature_march': None,
                  'temperature_april': None,
                  'temperature_may': None,
                  'temperature_june': None,
                  'temperature_july': None,
                  'temperature_august': None,
                  'temperature_september': None,
                  'temperature_october': None,
                  'temperature_november': None,
                  'temperature_december': None,
                  'outlet_temperature': 80,
                  'integration_scheme_name': 'NP_IE_ACS',
                  'integration_scheme_initials': 'NP_IE_ACS',
                  'aperture_area': 2,
                  'coll_n0': 0.75,
                  'coll_a1': 4,
                  'coll_a2': 0.003,
                  'coll_iam': 0.2,
                  'coll_test_flow': 0.07,
                  'coll_rows': 3,
                  'colls_per_row': 5,
                  'coll_tilt': 35,
                  'coll_azimuth': 0,
                  'field_mass_flow': 20,
                  'field_mass_flow_units': 'L/min',
                  'fluid': 'Agua',
                  'tank_volume': 5,
                  'tank_AR': 2.5,
                  'tank_material': None,
                  'tank_insulation_material': None,
                  'HX_eff': 0.8,
                  'fuel_cost': 1300,
                  'fuel_cost_units': '$/L'}


Result = simulate_system(system_params)

Heat_Map(Result['demanded_water_temperature'], 'Temperatura del agua demandada [°C]', 'Demanded_Temp')
Heat_Map(Result['useful_power_by_solar_field'], 'Calor aportado por campo solar [W]', 'Solar_Heat')
Heat_Map(Result['solar_frac'], 'Fracción solar', 'Solar_Frac')