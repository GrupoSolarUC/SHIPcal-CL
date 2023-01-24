# -*- coding: utf-8 -*-
"""
Created on Sun Nov  6 10:35:24 2022

@author: adria
"""

from Class_Solar_Field import Solar_Field
from Class_Heat_Exchanger import Heat_Exchanger
from Class_Storage import Storage_Tank
from Class_Properties import Properties
from simulation_functions import zenith_function, azimuth_function, monthly_flows, extract_weather_data, convert_power, convert_flow, convert_pressure, convert_demand, compute_cost_per_kJ, time_from_string, week_day_from_day, month_from_day
from math import sin, cos, pi
from pvlib.solarposition import get_solarposition
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d


def simulate_system(system_params, time_step = 0.1, tolerance = 1e-6, max_iterations = 100):
    
    
    latitude = -33.4981
    longitude = -70.6077
    weather_data = pd.read_csv('Weather_Data.csv')
    year_list, DNI, GHI, DHI, temp = extract_weather_data(weather_data)
    
    albedo = 0.25
    sky_diff = 0.5*( 1 + cos(system_params['coll_tilt']*pi/180) )*np.array(DHI)
    ground_diff = 0.5*albedo*( 1 - cos(system_params['coll_tilt']*pi/180) )*np.array(GHI)
    
    time_january = pd.date_range(str(year_list[0]) + '-01-01', str(year_list[0]) + '-02-01', inclusive='left', freq='min', tz = 'Etc/GMT+3')
    time_february = pd.date_range(str(year_list[1]) + '-02-01', str(year_list[1]) + '-03-01', inclusive='left', freq='min', tz = 'Etc/GMT+3')
    time_march = pd.date_range(str(year_list[2]) + '-03-01', str(year_list[2]) + '-04-01', inclusive='left', freq='min', tz = 'Etc/GMT+3')
    time_april = pd.date_range(str(year_list[3]) + '-04-01', str(year_list[3]) + '-05-01', inclusive='left', freq='min', tz = 'Etc/GMT+3')
    time_may = pd.date_range(str(year_list[4]) + '-05-01', str(year_list[4]) + '-06-01', inclusive='left', freq='min', tz = 'Etc/GMT+3')
    time_june = pd.date_range(str(year_list[5]) + '-06-01', str(year_list[5]) + '-07-01', inclusive='left', freq='min', tz = 'Etc/GMT+3')
    time_july = pd.date_range(str(year_list[6]) + '-07-01', str(year_list[6]) + '-08-01', inclusive='left', freq='min', tz = 'Etc/GMT+3')
    time_august = pd.date_range(str(year_list[7]) + '-08-01', str(year_list[7]) + '-09-01', inclusive='left', freq='min', tz = 'Etc/GMT+3')
    time_september = pd.date_range(str(year_list[8]) + '-09-01', str(year_list[8]) + '-10-01', inclusive='left', freq='min', tz = 'Etc/GMT+3')
    time_october = pd.date_range(str(year_list[9]) + '-10-01', str(year_list[9]) + '-11-01', inclusive='left', freq='min', tz = 'Etc/GMT+3')
    time_november = pd.date_range(str(year_list[10]) + '-11-01', str(year_list[10]) + '-12-01', inclusive='left', freq='min', tz = 'Etc/GMT+3')
    time_december = pd.date_range(str(year_list[11]) + '-12-01', str(year_list[11] + 1) + '-01-01', inclusive='left', freq='min', tz = 'Etc/GMT+3')

    if len(time_february) == 41760:
        time_february = pd.date_range(str(year_list[1]) + '-02-01', str(year_list[1]) + '-02-29', inclusive='left', freq='min', tz = 'Etc/GMT+3')
    
    solpos_january = get_solarposition(time_january, latitude, longitude)
    azimuth_january = solpos_january['azimuth']
    zenith_january = solpos_january['zenith']
    
    solpos_february = get_solarposition(time_february, latitude, longitude)
    azimuth_february = solpos_february['azimuth']
    zenith_february = solpos_february['zenith']
    
    solpos_march = get_solarposition(time_march, latitude, longitude)
    azimuth_march = solpos_march['azimuth']
    zenith_march = solpos_march['zenith']
    
    solpos_april = get_solarposition(time_april, latitude, longitude)
    azimuth_april = solpos_april['azimuth']
    zenith_april = solpos_april['zenith']
    
    solpos_may = get_solarposition(time_may, latitude, longitude)
    azimuth_may = solpos_may['azimuth']
    zenith_may = solpos_may['zenith']
    
    solpos_june = get_solarposition(time_june, latitude, longitude)
    azimuth_june = solpos_june['azimuth']
    zenith_june = solpos_june['zenith']
    
    solpos_july = get_solarposition(time_july, latitude, longitude)
    azimuth_july = solpos_july['azimuth']
    zenith_july = solpos_july['zenith']
    
    solpos_august = get_solarposition(time_august, latitude, longitude)
    azimuth_august = solpos_august['azimuth']
    zenith_august = solpos_august['zenith']
    
    solpos_september = get_solarposition(time_september, latitude, longitude)
    azimuth_september = solpos_september['azimuth']
    zenith_september = solpos_september['zenith']
    
    solpos_october = get_solarposition(time_october, latitude, longitude)
    azimuth_october = solpos_october['azimuth']
    zenith_october = solpos_october['zenith']
    
    solpos_november = get_solarposition(time_november, latitude, longitude)
    azimuth_november = solpos_november['azimuth']
    zenith_november = solpos_november['zenith']
    
    solpos_december = get_solarposition(time_december, latitude, longitude)
    azimuth_december = solpos_december['azimuth']
    zenith_december = solpos_december['zenith']
    
    solar_azimuth_list = (azimuth_january.tolist() +
                          azimuth_february.tolist() +
                          azimuth_march.tolist() +
                          azimuth_april.tolist() +
                          azimuth_may.tolist() +
                          azimuth_june.tolist() +
                          azimuth_july.tolist() +
                          azimuth_august.tolist() +
                          azimuth_september.tolist() +
                          azimuth_october.tolist() +
                          azimuth_november.tolist() +
                          azimuth_december.tolist() )
    
    solar_zenith_list = (zenith_january.tolist() +
                         zenith_february.tolist() +
                         zenith_march.tolist() +
                         zenith_april.tolist() +
                         zenith_may.tolist() +
                         zenith_june.tolist() +
                         zenith_july.tolist() +
                         zenith_august.tolist() +
                         zenith_september.tolist() +
                         zenith_october.tolist() +
                         zenith_november.tolist() +
                         zenith_december.tolist() )
    
    compute_azimuth = azimuth_function(solar_azimuth_list)
    compute_zenith = zenith_function(solar_zenith_list)
    
    time = range(8761)

    sky_diff_func = interp1d(time, sky_diff)
    ground_diff_func = interp1d(time, ground_diff)
    DNI_func = interp1d(time, DNI)
    T_amb_func = interp1d(time, temp)
    
    
    temp_january = np.mean(temp[:744])
    temp_february = np.mean(temp[744:1416])
    temp_march = np.mean(temp[1416:2160])
    temp_april = np.mean(temp[2160:2880])
    temp_may = np.mean(temp[2880:3624])
    temp_june = np.mean(temp[3624:4344])
    temp_july = np.mean(temp[4344:5088])
    temp_august = np.mean(temp[5088:5832])
    temp_september = np.mean(temp[5832:6552])
    temp_october = np.mean(temp[6552:7296])
    temp_november = np.mean(temp[7296:8016])
    temp_december = np.mean(temp[8016:8760])
    T_month_list = [temp_january,
                    temp_february,
                    temp_march,
                    temp_april,
                    temp_may,
                    temp_june,
                    temp_july,
                    temp_august,
                    temp_september,
                    temp_october,
                    temp_november,
                    temp_december ]
    T_amb_ann = np.mean(temp)
    delta_T_amb = ( max(T_month_list) - min(T_month_list) )/2
    delta_T_offset = 3.3333333
    T_ref = 6.6666667
    K1 = 0.4
    K2 = 0.018
    K3 = 35*pi/180
    K4 = -3.1416e-4
    delta_T_mains = (K1 + K2*(T_amb_ann - T_ref))*delta_T_amb
    phi_lag = K3 + K4*(T_amb_ann - T_ref)
    phi_amb = (104.8 + 180)*pi/180
    T_mains_avg = T_amb_ann + delta_T_offset
    def T_mains_func(t):
        '''
        Función que recibe el instante del año (en horas) y retorna la temperatura estimada para el agua de la red.

        Parámetros:
            - t: Instante del año (en horas)

        Retorna:
            - Temperatura estimada para el agua de la red (en °C)
        '''
        return T_mains_avg + delta_T_mains*sin(2*pi*t/8760 - phi_lag - phi_amb)
    
    monthly_mains_temperatures = {'January': np.mean( [ T_mains_func(t) for t in range(744) ] ),
                                  'February': np.mean( [ T_mains_func(t) for t in range(744, 1416) ] ),
                                  'March': np.mean( [ T_mains_func(t) for t in range(1416, 2160) ] ),
                                  'April': np.mean( [ T_mains_func(t) for t in range(2160, 2880) ] ),
                                  'May': np.mean( [ T_mains_func(t) for t in range(2880, 3624) ] ),
                                  'June': np.mean( [ T_mains_func(t) for t in range(3624, 4344) ] ),
                                  'July': np.mean( [ T_mains_func(t) for t in range(4344, 5088) ] ),
                                  'August': np.mean( [ T_mains_func(t) for t in range(5088, 5832) ] ),
                                  'September': np.mean( [ T_mains_func(t) for t in range(5832, 6552) ] ),
                                  'October': np.mean( [ T_mains_func(t) for t in range(6552, 7296) ] ),
                                  'November': np.mean( [ T_mains_func(t) for t in range(7296, 8016) ] ),
                                  'December': np.mean( [ T_mains_func(t) for t in range(8016, 8760) ] ) }
    
    if system_params['closed_system']:
        if system_params['monthly_return_temperature']:
            monthly_inlet_temperatures = {'January': system_params['temperature_january'],
                                          'February': system_params['temperature_february'],
                                          'March': system_params['temperature_march'],
                                          'April': system_params['temperature_april'],
                                          'May': system_params['temperature_may'],
                                          'June': system_params['temperature_june'],
                                          'July': system_params['temperature_july'],
                                          'August': system_params['temperature_august'],
                                          'September': system_params['temperature_september'],
                                          'October': system_params['temperature_october'],
                                          'November': system_params['temperature_november'],
                                          'December': system_params['temperature_december'] }
        else:
            monthly_inlet_temperatures = {'January': system_params['return_inlet_temperature'],
                                          'February': system_params['return_inlet_temperature'],
                                          'March': system_params['return_inlet_temperature'],
                                          'April': system_params['return_inlet_temperature'],
                                          'May': system_params['return_inlet_temperature'],
                                          'June': system_params['return_inlet_temperature'],
                                          'July': system_params['return_inlet_temperature'],
                                          'August': system_params['return_inlet_temperature'],
                                          'September': system_params['return_inlet_temperature'],
                                          'October': system_params['return_inlet_temperature'],
                                          'November': system_params['return_inlet_temperature'],
                                          'December': system_params['return_inlet_temperature'] }
    else:
        monthly_inlet_temperatures = monthly_mains_temperatures
    
    field_fluid = system_params['fluid']
    
    if field_fluid == 'Agua':
        PropsField = Properties(0, 5)
    else:
        glycol_percentage = int(system_params['fluid'][7:9])
        PropsField = Properties(glycol_percentage, 5)
    
    PropsHeatedFluid = Properties(0, 2)
    PropsBoiler = Properties(0, convert_pressure(system_params['boiler_pressure'], system_params['boiler_pressure_units']))
    
    work_day_list = []
    if system_params['demand_monday']:
        work_day_list.append('Monday')
    if system_params['demand_tuesday']:
        work_day_list.append('Tuesday')
    if system_params['demand_wednesday']:
        work_day_list.append('Wednesday')
    if system_params['demand_thursday']:
        work_day_list.append('Thursday')
    if system_params['demand_friday']:
        work_day_list.append('Friday')
    if system_params['demand_saturday']:
        work_day_list.append('Saturday')
    if system_params['demand_sunday']:
        work_day_list.append('Sunday')
    
    monthly_demands = {'January': convert_demand(system_params['demand_january'], system_params['yearly_demand_unit'], system_params['fuel_name'], system_params['boiler_type'], system_params['boiler_efficiency']),
                       'February': convert_demand(system_params['demand_february'], system_params['yearly_demand_unit'], system_params['fuel_name'], system_params['boiler_type'], system_params['boiler_efficiency']),
                       'March': convert_demand(system_params['demand_march'], system_params['yearly_demand_unit'], system_params['fuel_name'], system_params['boiler_type'], system_params['boiler_efficiency']),
                       'April': convert_demand(system_params['demand_april'], system_params['yearly_demand_unit'], system_params['fuel_name'], system_params['boiler_type'], system_params['boiler_efficiency']),
                       'May': convert_demand(system_params['demand_may'], system_params['yearly_demand_unit'], system_params['fuel_name'], system_params['boiler_type'], system_params['boiler_efficiency']),
                       'June': convert_demand(system_params['demand_june'], system_params['yearly_demand_unit'], system_params['fuel_name'], system_params['boiler_type'], system_params['boiler_efficiency']),
                       'July': convert_demand(system_params['demand_july'], system_params['yearly_demand_unit'], system_params['fuel_name'], system_params['boiler_type'], system_params['boiler_efficiency']),
                       'August': convert_demand(system_params['demand_august'], system_params['yearly_demand_unit'], system_params['fuel_name'], system_params['boiler_type'], system_params['boiler_efficiency']),
                       'September': convert_demand(system_params['demand_september'], system_params['yearly_demand_unit'], system_params['fuel_name'], system_params['boiler_type'], system_params['boiler_efficiency']),
                       'October': convert_demand(system_params['demand_october'], system_params['yearly_demand_unit'], system_params['fuel_name'], system_params['boiler_type'], system_params['boiler_efficiency']),
                       'November': convert_demand(system_params['demand_november'], system_params['yearly_demand_unit'], system_params['fuel_name'], system_params['boiler_type'], system_params['boiler_efficiency']),
                       'December': convert_demand(system_params['demand_december'], system_params['yearly_demand_unit'], system_params['fuel_name'], system_params['boiler_type'], system_params['boiler_efficiency'])}
    
    t_start = time_from_string(system_params['operation_start'])
    t_end = time_from_string(system_params['operation_end'])
    
    flows = monthly_flows(monthly_demands, monthly_inlet_temperatures, work_day_list, t_start, t_end, system_params['outlet_temperature'], PropsHeatedFluid)
    Field_mass_flow = convert_flow(system_params['field_mass_flow'], system_params['field_mass_flow_units'], PropsField.glycol_percentage)
    
    if latitude < 49:
        def corrected_time(t):
            '''
            Función que recibe la hora del año (entre 0 y 8760) y retorna la hora considerando el cambio de hora (entre 0 y 8760).
            
            Esta función considera que su argumento es la hora según el horario UTC-3, y en caso de encontrarse en el rango de tiempo durante el cual la hora oficial es UTC-4, resta 1 al valor recibido.
            
            En caso de que la simulación se realice en la región de Magallanes (el criterio para estar en la región de Magallanes es que la latitud sea más al sur que -49°), la función SIEMPRE retorna el mismo valor recibido como argumento.
            
            Se considera que los cambios de hora se realizan los días 7 de abril y 8 de septiembre a las 24:00 horas.
    
            Prámetros:
                - t: Hora del año de acuerdo al horario UTC-3
    
            Retorna:
                - Hora del año de acuerdo al horario UTC-3 o UTC-4, según corresponda.
            '''
            if t >= 24*97 and t - 1 < 24*251:
                return t - 1
            else:
                return t
    else:
        def corrected_time(t):
            '''
            Función que recibe la hora del año (entre 0 y 8760) y retorna la hora considerando el cambio de hora (entre 0 y 8760).
            
            Esta función considera que su argumento es la hora según el horario UTC-3, y en caso de encontrarse en el rango de tiempo durante el cual la hora oficial es UTC-4, resta 1 al valor recibido.
            
            En caso de que la simulación se realice en la región de Magallanes (el criterio para estar en la región de Magallanes es que la latitud sea más al sur que -49°), la función SIEMPRE retorna el mismo valor recibido como argumento.
            
            Se considera que los cambios de hora se realizan los días 7 de abril y 8 de septiembre a las 24:00 horas.
    
            Prámetros:
                - t: Hora del año de acuerdo al horario UTC-3
    
            Retorna:
                - Hora del año de acuerdo al horario UTC-3 o UTC-4, según corresponda.
            '''
            return t
    
    cost_per_kJ = compute_cost_per_kJ(system_params['fuel_name'], system_params['fuel_cost'], system_params['fuel_cost_units'], system_params['boiler_type'])
    
    scheme_list_1 = ['NP_IE_ACS', 'NS_L_SI', 'NS_L_PI', 'SAM']
    scheme_list_2 = ['NS_L_PD', 'NS_L_CA_MU', 'NS_L_CA_1', 'NS_L_CA_2']
    
    ###########################################################################################################
    ###########################################################################################################
    
    if system_params['integration_scheme_initials'] in scheme_list_1 or system_params['integration_scheme_name'] in scheme_list_1:
        
        Tank = Storage_Tank(system_params['tank_volume'], system_params['tank_AR'], 5, 5, 5,
                            10, None, None, 10, 1,
                            system_params['HX_eff'], system_params['HX_eff'], 7, 9, 2, 4,
                            PropsHeatedFluid, PropsField, PropsBoiler, {i: 20 for i in range(1,11)} )
        
        Field = Solar_Field(system_params['aperture_area'], system_params['coll_n0'],
                            system_params['coll_a1'], system_params['coll_a2'],
                            system_params['coll_iam'], system_params['coll_test_flow']*3600, Properties(0, 2),
                            system_params['coll_tilt'], system_params['coll_azimuth'],
                            system_params['coll_rows'], system_params['colls_per_row'], PropsField)
        
        boiler_nominal_power = convert_power(system_params['boiler_nominal_power'], system_params['boiler_nominal_power_units'])
        
        rad_min = (system_params['coll_a1']*15 + system_params['coll_a2']*15**2)/system_params['coll_n0']
        
        h_in_Boiler_HX = PropsBoiler.T_to_h(system_params['outlet_temperature'])
        flow_boiler = 2*max([ flows[key] for key in flows ])
        
        t = 8760 - 24*10
        while t < 8760:
            day = int(t/24) + 1
            month = month_from_day(day)
            week_day = week_day_from_day(day)
            day_time = corrected_time(t)%24
            
            sky_diffuse = float(sky_diff_func(t))
            ground_diffuse = float(ground_diff_func(t))
            DNIrr = float(DNI_func(t))
            T_amb = float(T_amb_func(t))
            T_mains = monthly_mains_temperatures[month]
            h_mains = PropsHeatedFluid.T_to_h(T_mains)
            h_mains_field = PropsField.T_to_h(T_mains)
            demanded_flow = flows[month]
            aoi, longi, trans = Field.incidence_angles(float(compute_zenith(t)), float(compute_azimuth(t)))
        
            if t_start == t_end and demanded_flow > 0 and week_day in work_day_list:
                system_operation = True
            elif t_start < t_end and day_time >= t_start and day_time <= t_end and demanded_flow > 0 and week_day in work_day_list:
                system_operation = True
            elif t_start > t_end and (day_time >= t_start or day_time <= t_end) and demanded_flow > 0 and week_day in work_day_list:
                system_operation = True
            else:
                system_operation = False
                solar_field_operation = False
            
            if system_operation:
                if aoi == None:
                    rad = sky_diffuse + ground_diffuse
                    if rad > 0:
                        IAM_eff = (sky_diffuse*Field.sky_diffuse_iam + ground_diffuse*Field.ground_diffuse_iam)/rad
                else:
                    beam_rad = DNIrr*cos(aoi)
                    beam_iam = Field.IAM_func(aoi, longi, trans)
                    rad = beam_rad + sky_diffuse + ground_diffuse
                    if rad > 0:
                        IAM_eff = (beam_rad*beam_iam + sky_diffuse*Field.sky_diffuse_iam + ground_diffuse*Field.ground_diffuse_iam)/rad
                if rad >= rad_min:
                    solar_field_operation = True
                else:
                    solar_field_operation = False
                if solar_field_operation:
                    h_in_solar_field = PropsField.T_to_h(40)
                    it = 0
                    while True:
                        Field_outputs = Field.compute_outputs(Field_mass_flow, h_in_solar_field, h_mains_field, rad, IAM_eff, T_amb)
                        Tank_outputs = Tank.compute_outputs(0, h_mains, demanded_flow, h_mains,
                                                            Field_mass_flow, Field_outputs['h_out'], flow_boiler, h_in_Boiler_HX,
                                                            T_amb, time_step)
                        it = it + 1
                        if abs(Tank_outputs['HX1_outlet_h'] - h_in_solar_field)/h_in_solar_field < tolerance:
                            maxIt = False
                            break
                        if it == 100:
                            maxIt = True
                            break
                        h_in_solar_field = Tank_outputs['HX1_outlet_h']
                else:
                    Tank_outputs = Tank.compute_outputs(0, h_mains, demanded_flow, h_mains,
                                                        None, None, flow_boiler, h_in_Boiler_HX,
                                                        T_amb, time_step)
            else:
                Tank_outputs = Tank.compute_outputs(0, h_mains, 0, h_mains,
                                                    None, None, None, None,
                                                    T_amb, time_step)
            Tank.update_temperature()
            t = np.round(t + time_step, 2)
        
        Result = []
        t = 0
        while t < 8760:
            day = int(t/24) + 1
            month = month_from_day(day)
            week_day = week_day_from_day(day)
            day_time = corrected_time(t)%24
            
            sky_diffuse = float(sky_diff_func(t))
            ground_diffuse = float(ground_diff_func(t))
            DNIrr = float(DNI_func(t))
            T_amb = float(T_amb_func(t))
            T_mains = monthly_mains_temperatures[month]
            h_mains = PropsHeatedFluid.T_to_h(T_mains)
            h_mains_field = PropsField.T_to_h(T_mains)
            demanded_flow = flows[month]
            aoi, longi, trans = Field.incidence_angles(float(compute_zenith(t)), float(compute_azimuth(t)))
        
            if t_start == t_end and demanded_flow > 0 and week_day in work_day_list:
                system_operation = True
            elif t_start < t_end and day_time >= t_start and day_time <= t_end and demanded_flow > 0 and week_day in work_day_list:
                system_operation = True
            elif t_start > t_end and (day_time >= t_start or day_time <= t_end) and demanded_flow > 0 and week_day in work_day_list:
                system_operation = True
            else:
                system_operation = False
                solar_field_operation = False
            
            if system_operation:
                if aoi == None:
                    rad = sky_diffuse + ground_diffuse
                    if rad > 0:
                        IAM_eff = (sky_diffuse*Field.sky_diffuse_iam + ground_diffuse*Field.ground_diffuse_iam)/rad
                else:
                    beam_rad = DNIrr*cos(aoi)
                    beam_iam = Field.IAM_func(aoi, longi, trans)
                    rad = beam_rad + sky_diffuse + ground_diffuse
                    if rad > 0:
                        IAM_eff = (beam_rad*beam_iam + sky_diffuse*Field.sky_diffuse_iam + ground_diffuse*Field.ground_diffuse_iam)/rad
                if rad >= rad_min:
                    solar_field_operation = True
                else:
                    solar_field_operation = False
                if solar_field_operation:
                    h_in_solar_field = PropsField.T_to_h(40)
                    it = 0
                    while True:
                        Field_outputs = Field.compute_outputs(Field_mass_flow, h_in_solar_field, h_mains_field, rad, IAM_eff, T_amb)
                        Tank_outputs = Tank.compute_outputs(0, h_mains, demanded_flow, h_mains,
                                                            Field_mass_flow, Field_outputs['h_out'], flow_boiler, h_in_Boiler_HX,
                                                            T_amb, time_step)
                        it = it + 1
                        if abs(Tank_outputs['HX1_outlet_h'] - h_in_solar_field)/h_in_solar_field < tolerance:
                            maxIt = False
                            break
                        if it == 100:
                            maxIt = True
                            break
                        h_in_solar_field = Tank_outputs['HX1_outlet_h']
                else:
                    Tank_outputs = Tank.compute_outputs(0, h_mains, demanded_flow, h_mains,
                                                        None, None, flow_boiler, h_in_Boiler_HX,
                                                        T_amb, time_step)
            else:
                Tank_outputs = Tank.compute_outputs(0, h_mains, 0, h_mains,
                                                    None, None, None, None,
                                                    T_amb, time_step)
            if solar_field_operation:
                wasted_heat = Field_outputs['Q_waste'] + Tank_outputs['Q_waste']
            else:
                wasted_heat = 0
                
            if system_operation:
                heat_rate_to_demand = demanded_flow*(Tank_outputs['outlet_2_h'] - h_mains)
            else:
                heat_rate_to_demand = 0
            
            Result.append({'t': t,
                           'total_irradiance': rad/3.6,
                           'system_operation': system_operation,
                           'solar_field_operation': solar_field_operation,
                           'T_mains': T_mains,
                           'demanded_water_temperature': PropsHeatedFluid.h_to_T(Tank_outputs['outlet_2_h']),
                           'total_power_delivered_to_demand': heat_rate_to_demand/3.6,
                           'useful_power_by_solar_field': Tank_outputs['HX1_Q']/3.6,
                           'wasted_heat': wasted_heat/3.6,
                           'useful_power_by_boiler': Tank_outputs['HX2_Q']/3.6,
                           'maximum_iterations_reached': maxIt,
                           'saturation_inside_tank': Tank_outputs['T_sat_reached']})
            Tank.update_temperature()
            t = np.round(t + time_step, 2)
        
        Result = { key: [ Result[i][key] for i in range(len(Result)) ] for key in Result[0] }
        
        solar_frac = []
        for i in range(len(Result['useful_power_by_boiler'])):
            if Result['useful_power_by_boiler'][i] > 0 or Result['useful_power_by_solar_field'][i] > 0:
                solar_frac.append( min ( [ Result['useful_power_by_solar_field'][i]/(Result['useful_power_by_boiler'][i] + Result['useful_power_by_solar_field'][i]), 1]) )
            else:
                solar_frac.append(np.nan)
                
        Result['solar_frac'] = solar_frac
        
    ###########################################################################################################
    ###########################################################################################################
            
    if system_params['integration_scheme_initials'] in scheme_list_2 or system_params['integration_scheme_name'] in scheme_list_2:
        
        Tank = Storage_Tank(system_params['tank_volume'], system_params['tank_AR'], 5, 5, 5,
                            10, 7, 10, 10, 1,
                            None, system_params['HX_eff'], None, None, 2, 4,
                            PropsHeatedFluid, None, PropsBoiler, {i: 20 for i in range(1,11)} )
        
        Field = Solar_Field(system_params['aperture_area'], system_params['coll_n0'],
                            system_params['coll_a1'], system_params['coll_a2'],
                            system_params['coll_iam'], system_params['coll_test_flow']*3600, Properties(0, 2),
                            system_params['coll_tilt'], system_params['coll_azimuth'],
                            system_params['coll_rows'], system_params['colls_per_row'], PropsField)
        
        HX = Heat_Exchanger(system_params['HX_eff'], PropsField, PropsHeatedFluid)
        
        boiler_nominal_power = convert_power(system_params['boiler_nominal_power'], system_params['boiler_nominal_power_units'])
        
        rad_min = (system_params['coll_a1']*15 + system_params['coll_a2']*15**2)/system_params['coll_n0']
        
        h_in_Boiler_HX = PropsBoiler.T_to_h(system_params['outlet_temperature'])
        flow_boiler = 2*max([ flows[key] for key in flows ])
        
        t = 8760 - 24*10
        while t < 8760:
            day = int(t/24) + 1
            month = month_from_day(day)
            week_day = week_day_from_day(day)
            day_time = corrected_time(t)%24
            
            sky_diffuse = float(sky_diff_func(t))
            ground_diffuse = float(ground_diff_func(t))
            DNIrr = float(DNI_func(t))
            T_amb = float(T_amb_func(t))
            T_mains = monthly_mains_temperatures[month]
            h_mains = PropsHeatedFluid.T_to_h(T_mains)
            h_mains_field = PropsField.T_to_h(T_mains)
            demanded_flow = flows[month]
            aoi, longi, trans = Field.incidence_angles(float(compute_zenith(t)), float(compute_azimuth(t)))
        
            if t_start == t_end and demanded_flow > 0 and week_day in work_day_list:
                system_operation = True
            elif t_start < t_end and day_time >= t_start and day_time <= t_end and demanded_flow > 0 and week_day in work_day_list:
                system_operation = True
            elif t_start > t_end and (day_time >= t_start or day_time <= t_end) and demanded_flow > 0 and week_day in work_day_list:
                system_operation = True
            else:
                system_operation = False
                solar_field_operation = False
            
            if system_operation:
                if aoi == None:
                    rad = sky_diffuse + ground_diffuse
                    if rad > 0:
                        IAM_eff = (sky_diffuse*Field.sky_diffuse_iam + ground_diffuse*Field.ground_diffuse_iam)/rad
                else:
                    beam_rad = DNIrr*cos(aoi)
                    beam_iam = Field.IAM_func(aoi, longi, trans)
                    rad = beam_rad + sky_diffuse + ground_diffuse
                    if rad > 0:
                        IAM_eff = (beam_rad*beam_iam + sky_diffuse*Field.sky_diffuse_iam + ground_diffuse*Field.ground_diffuse_iam)/rad
                if rad >= rad_min:
                    solar_field_operation = True
                else:
                    solar_field_operation = False
                if solar_field_operation:
                    previous_enthalpies = {'h_in_solar_field': PropsField.T_to_h(40), 'h_in_solar_HX_load': PropsHeatedFluid.T_to_h(40) }
                    it = 0
                    while True:
                        Field_outputs = Field.compute_outputs(Field_mass_flow, previous_enthalpies['h_in_solar_field'], h_mains_field, rad, IAM_eff, T_amb)
                        HX_outputs = HX.compute_outputs(Field_mass_flow, Field_outputs['h_out'], Field_mass_flow, previous_enthalpies['h_in_solar_HX_load'], h_mains)
                        Tank_outputs = Tank.compute_outputs(Field_mass_flow, HX_outputs['h_out_load'], demanded_flow, h_mains, None, None, flow_boiler, h_in_Boiler_HX, T_amb, time_step)
                        it = it + 1
                        if ( abs(HX_outputs['h_out_source'] - previous_enthalpies['h_in_solar_field'])/previous_enthalpies['h_in_solar_field'] < tolerance and
                             abs(Tank_outputs['outlet_1_h'] - previous_enthalpies['h_in_solar_HX_load'])/previous_enthalpies['h_in_solar_HX_load'] < tolerance):
                            maxIt = False
                            break
                        if it == max_iterations:
                            maxIt = True
                            break
                        previous_enthalpies = {'h_in_solar_field': HX_outputs['h_out_source'], 'h_in_solar_HX_load': Tank_outputs['outlet_1_h'] }
                else:
                    Tank_outputs = Tank.compute_outputs(0, h_mains, demanded_flow, h_mains, None, None, flow_boiler, h_in_Boiler_HX, T_amb, time_step)
            else:
                Tank_outputs = Tank.compute_outputs(0, h_mains, 0, h_mains, None, None, None, None, T_amb, time_step)
            Tank.update_temperature()
            t = np.round(t + time_step, 2)
        
        Result = []
        t = 0
        while t < 8760:
            day = int(t/24) + 1
            month = month_from_day(day)
            week_day = week_day_from_day(day)
            day_time = corrected_time(t)%24
            
            sky_diffuse = float(sky_diff_func(t))
            ground_diffuse = float(ground_diff_func(t))
            DNIrr = float(DNI_func(t))
            T_amb = float(T_amb_func(t))
            T_mains = monthly_mains_temperatures[month]
            h_mains = PropsHeatedFluid.T_to_h(T_mains)
            h_mains_field = PropsField.T_to_h(T_mains)
            demanded_flow = flows[month]
            aoi, longi, trans = Field.incidence_angles(float(compute_zenith(t)), float(compute_azimuth(t)))
        
            if t_start == t_end and demanded_flow > 0 and week_day in work_day_list:
                system_operation = True
            elif t_start < t_end and day_time >= t_start and day_time <= t_end and demanded_flow > 0 and week_day in work_day_list:
                system_operation = True
            elif t_start > t_end and (day_time >= t_start or day_time <= t_end) and demanded_flow > 0 and week_day in work_day_list:
                system_operation = True
            else:
                system_operation = False
                solar_field_operation = False
            
            if system_operation:
                if aoi == None:
                    rad = sky_diffuse + ground_diffuse
                    if rad > 0:
                        IAM_eff = (sky_diffuse*Field.sky_diffuse_iam + ground_diffuse*Field.ground_diffuse_iam)/rad
                else:
                    beam_rad = DNIrr*cos(aoi)
                    beam_iam = Field.IAM_func(aoi, longi, trans)
                    rad = beam_rad + sky_diffuse + ground_diffuse
                    if rad > 0:
                        IAM_eff = (beam_rad*beam_iam + sky_diffuse*Field.sky_diffuse_iam + ground_diffuse*Field.ground_diffuse_iam)/rad
                if rad >= rad_min:
                    solar_field_operation = True
                else:
                    solar_field_operation = False
                if solar_field_operation:
                    previous_enthalpies = {'h_in_solar_field': PropsField.T_to_h(40), 'h_in_solar_HX_load': PropsHeatedFluid.T_to_h(40) }
                    it = 0
                    while True:
                        Field_outputs = Field.compute_outputs(Field_mass_flow, previous_enthalpies['h_in_solar_field'], h_mains_field, rad, IAM_eff, T_amb)
                        HX_outputs = HX.compute_outputs(Field_mass_flow, Field_outputs['h_out'], Field_mass_flow, previous_enthalpies['h_in_solar_HX_load'], h_mains)
                        Tank_outputs = Tank.compute_outputs(Field_mass_flow, HX_outputs['h_out_load'], demanded_flow, h_mains, None, None, flow_boiler, h_in_Boiler_HX, T_amb, time_step)
                        it = it + 1
                        if ( abs(HX_outputs['h_out_source'] - previous_enthalpies['h_in_solar_field'])/previous_enthalpies['h_in_solar_field'] < tolerance and
                             abs(Tank_outputs['outlet_1_h'] - previous_enthalpies['h_in_solar_HX_load'])/previous_enthalpies['h_in_solar_HX_load'] < tolerance):
                            maxIt = False
                            break
                        if it == max_iterations:
                            maxIt = True
                            break
                        previous_enthalpies = {'h_in_solar_field': HX_outputs['h_out_source'], 'h_in_solar_HX_load': Tank_outputs['outlet_1_h'] }
                else:
                    Tank_outputs = Tank.compute_outputs(0, h_mains, demanded_flow, h_mains, None, None, flow_boiler, h_in_Boiler_HX, T_amb, time_step)
            else:
                Tank_outputs = Tank.compute_outputs(0, h_mains, 0, h_mains, None, None, None, None, T_amb, time_step)
            
            if solar_field_operation:
                wasted_heat = Field_outputs['Q_waste'] + HX_outputs['Q_waste']
            else:
                wasted_heat = 0
                
            if system_operation:
                heat_rate_to_demand = demanded_flow*(Tank_outputs['outlet_2_h'] - h_mains)
            else:
                heat_rate_to_demand = 0
            
            Result.append({'t': t,
                           'total_irradiance': rad/3.6,
                           'system_operation': system_operation,
                           'solar_field_operation': solar_field_operation,
                           'T_mains': T_mains,
                           'demanded_water_temperature': PropsHeatedFluid.h_to_T(Tank_outputs['outlet_2_h']),
                           'total_power_delivered_to_demand': heat_rate_to_demand/3.6,
                           'useful_power_by_solar_field': HX_outputs['Q_useful']/3.6,
                           'wasted_heat': wasted_heat/3.6,
                           'useful_power_by_boiler': Tank_outputs['HX2_Q']/3.6,
                           'maximum_iterations_reached': maxIt,
                           'Tank_T_sat_reached': Tank_outputs['T_sat_reached']})
                
            Tank.update_temperature()
            t = np.round(t + time_step, 2)
        
        Result = {key: [ Result[i][key] for i in range(len(Result)) ] for key in Result[0]}
        
        solar_frac = []
        for i in range(len(Result['useful_power_by_boiler'])):
            if Result['useful_power_by_boiler'][i] > 0 or Result['useful_power_by_solar_field'][i] > 0:
                solar_frac.append( min ( [ Result['useful_power_by_solar_field'][i]/(Result['useful_power_by_boiler'][i] + Result['useful_power_by_solar_field'][i]), 1]) )
            else:
                solar_frac.append(np.nan)
                
        Result['solar_frac'] = solar_frac
        
        
    return Result




















