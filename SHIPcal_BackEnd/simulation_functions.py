#!/usr/bin/env python
# coding: utf-8



import numpy as np
import pandas as pd
import pvlib as pv
from scipy import interpolate
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import interp1d
from scipy.optimize import minimize, Bounds

def zenith_function(zenith_list):
    assert len(zenith_list)%8760 == 0
    factor = int(len(zenith_list)/8760)
    time = [ t/factor for t in range(len(zenith_list)) ]
    time.append(8760)
    zenith_list.append(zenith_list[0])
    return interp1d(time, zenith_list)

def azimuth_function(azimuth_list):
    assert len(azimuth_list)%8760 == 0
    factor = int(len(azimuth_list)/8760)
    time = [ t/factor for t in range(len(azimuth_list)) ]
    i = 0
    while i <  len(azimuth_list) - 1:
        if azimuth_list[i + 1] - azimuth_list[i] > 200:
            if azimuth_list[i] == 0 and azimuth_list[i + 1] == 360:
                i = i + 1
            else:
                delta_t = time[i + 1] - time[i]
                delta_t_1 = azimuth_list[i]*delta_t/(azimuth_list[i] + (360 - azimuth_list[i + 1]))
                time.insert(i + 1, time[i] + delta_t_1)
                azimuth_list.insert(i + 1, 0)
                time.insert(i + 2, time[i] + delta_t_1)
                azimuth_list.insert(i + 2, 360)
                i = i + 3
        elif azimuth_list[i] - azimuth_list[i + 1] > 200:
            if azimuth_list[i] == 360 and azimuth_list[i + 1] == 0:
                i = i + 1
            else:
                delta_t = time[i + 1] - time[i]
                delta_t_1 = (360 - azimuth_list[i])*delta_t/((360 - azimuth_list[i]) + azimuth_list[i + 1])
                time.insert(i + 1, time[i] + delta_t_1)
                azimuth_list.insert(i + 1, 360)
                time.insert(i + 2, time[i] + delta_t_1)
                azimuth_list.insert(i + 2, 0)
                i = i + 3
        else:
            i = i + 1
    time.append(8760)
    azimuth_list.append(azimuth_list[0])
    return interp1d(time, azimuth_list)

def Heat_Map(Curve, title, filename):
    '''
    Función para generar un "mapa" de alguna variable (con las horas del día en el eje vertical y los días del año en el eje horizontal)

    Parámetros:
        - Curve : Lista de datos de la variable a plotear. Debe ser de un largo divisible en 365, de forma que para cada día del año haya la misma cantidad de datos.
        - title : Título de la figura
        - filename : Nombre del archivo (imagen) con el que se guardará el plot.

    Retorna:
        - None
    '''
    if len(Curve)%365 != 0:
        print('Número de elementos en el vector de resultados no es válido')
        return
    time_steps_per_day = int( len( Curve )/365 )
    plt.ioff()
    plt.rc('xtick', labelsize=14)
    plt.rc('ytick', labelsize=14)
    yticks = [ 0, int(time_steps_per_day/6), int(time_steps_per_day/3),
               int(time_steps_per_day/2), int(time_steps_per_day*2/3),
               int(time_steps_per_day*5/6), time_steps_per_day]
    ylabels = [0,4,8,12,16,20,24]
    matrix = np.reshape(Curve, ( 365, time_steps_per_day ) )
    matrix = np.transpose(matrix)
    fig = plt.figure(figsize = (8,12), dpi = 200)
    ax = plt.gca()
    cmap = plt.cm.get_cmap("jet").copy()
    cmap.set_bad('white',1.)
    im = ax.matshow(matrix, cmap = cmap, origin = 'lower')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.2)
    plt.colorbar(im, cax=cax)
    ax.set_yticks(yticks, ylabels)
    ax.xaxis.set_ticks_position('bottom')
    ax.xaxis.set_label_position('bottom')
    ax.set_ylabel('Hora del día', fontsize = 14)
    ax.set_xlabel('Día del año', fontsize = 14)
    ax.set_title(title, fontsize = 16, loc = 'left')
    plt.savefig(filename+'.jpg', bbox_inches='tight')

def monthly_work_hours(month, t_start, t_end, work_day_list):
    if t_start == t_end:
        daily_working_hours = 24
    else:
        if t_end > t_start:
            daily_working_hours = t_end - t_start
        else:
            daily_working_hours = 24 - t_start + t_end
    work_hours = 0
    for day in range(1,366):
        if month_from_day(day) == month and week_day_from_day(day) in work_day_list:
            work_hours = work_hours + daily_working_hours
    return work_hours

def monthly_flows(monthly_demands, monthly_temperatures, work_day_list, t_start, t_end, T_set, fluid_properties):
    flows = {}
    for month in monthly_demands:
        delta_h = fluid_properties.T_to_h(T_set) - fluid_properties.T_to_h(monthly_temperatures[month])
        if monthly_work_hours(month, t_start, t_end, work_day_list) == 0:
            flows[month] = 0
        else:
            m_month = monthly_demands[month]/delta_h
            flows[month] = m_month/monthly_work_hours(month, t_start, t_end, work_day_list)
    return flows

def extract_weather_data(dataframe):
    years = dataframe['Year'].values
    year_january = years[0]
    year_february = years[744]
    year_march = years[1416]
    year_april = years[2160]
    year_may = years[2880]
    year_june = years[3624]
    year_july = years[4344]
    year_august = years[5088]
    year_september = years[5832]
    year_october = years[6552]
    year_november = years[7296]
    year_december = years[8016]
    DNI = dataframe['DNI'].values
    GHI = dataframe['GHI'].values
    DHI = dataframe['DHI'].values
    temp = dataframe['Tdry']
    DNI = [ DNI[8759] ] + DNI[:8759].tolist()
    GHI = [ GHI[8759] ] + GHI[:8759].tolist()
    DHI = [ DHI[8759] ] + DHI[:8759].tolist() 
    temp = [ temp[8759] ] + temp[:8759].tolist()
    DNI.append(DNI[0])
    GHI.append(GHI[0])
    DHI.append(DHI[0])
    temp.append(temp[0])
    assert len(DNI) == 8761 and len(GHI) == 8761 and len(DHI) == 8761 and len(temp) == 8761
    DNI = [ 3.6*DNI[i] for i in range(len(DNI)) ]
    GHI = [ 3.6*GHI[i] for i in range(len(GHI)) ]
    DHI = [ 3.6*DHI[i] for i in range(len(DHI)) ]
    year_list = [year_january,
                 year_february,
                 year_march,
                 year_april,
                 year_may,
                 year_june,
                 year_july,
                 year_august,
                 year_september,
                 year_october,
                 year_november,
                 year_december ]
    return year_list, DNI, GHI, DHI, temp

def convert_power(original_value, original_units):
    if original_units == 'kJ/h':
        return original_value
    if original_units == 'kW':
        return original_value*3600
    if original_units == 'Ton/h':
        return original_value*2257000

def convert_flow(original_value, original_units, glycol_percentage):
    if original_units == 'kg/h':
        return original_value
    if original_units == 'kg/min':
        return original_value*60
    if original_units == 'kg/s':
        return original_value*3600
    glycol_fraction = glycol_percentage/100
    density = 1.11/( (1 - glycol_fraction)*1.11 + glycol_fraction )
    if original_units == 'L/h':
        return density*original_value
    if original_units == 'L/min':
        return density*original_value*60
    if original_units == 'L/s':
        return density*original_value*3600
    if original_units == 'gpm':
        return 227.125*density*original_value
    
def convert_pressure(original_value, original_units):
    if original_units == 'bar':
        return original_value
    if original_units == 'MPa':
        return original_value*10
    if original_units == 'psi':
        return original_value*0.06895

def convert_demand(original_value, original_units, fuel, boiler_type, boiler_efficiency):
    if original_units == 'MJ':
        return original_value*1000*boiler_efficiency
    if original_units == 'kWh':
        return original_value*3600*boiler_efficiency
    if original_units == 'MWh':
        return original_value*3600000*boiler_efficiency
    if original_units == 'BTU':
        return original_value*1.055*boiler_efficiency
    if original_units == 'kcal':
        return original_value*4.184*boiler_efficiency
    if fuel == 'Electricidad':
        raise ValueError('El consumo eléctrico solo puede expresarse en unidades de energía (MJ, kWh, MWh, BTU o kcal); no en kg, L ni m3')
    if fuel == 'Propano':
        fuel_density = 1.83
        if boiler_type == 'Condensación líquido' or boiler_type == 'Condensación vapor':
            fuel_heat_value = 50426
        else:
            fuel_heat_value = 46367
    if fuel == 'Butano':
        fuel_density = 2.52
        if boiler_type == 'Condensación líquido' or boiler_type == 'Condensación vapor':
            fuel_heat_value = 49652
        else:
            fuel_heat_value = 45765
    if fuel == 'Biomasa':
        fuel_density = 800
        if boiler_type == 'Condensación líquido' or boiler_type == 'Condensación vapor':
            fuel_heat_value = 19000
        else:
            fuel_heat_value = 17000
    if fuel == 'Diesel':
        fuel_density = 850
        if boiler_type == 'Condensación líquido' or boiler_type == 'Condensación vapor':
            fuel_heat_value = 45600
        else:
            fuel_heat_value = 42600
    if original_units == 'L':
        kJ_per_litre = fuel_density*fuel_heat_value/1000
        return kJ_per_litre*original_value*boiler_efficiency
    if original_units == 'm3':
        kJ_per_m3 = fuel_density*fuel_heat_value
        return kJ_per_m3*original_value*boiler_efficiency
    if original_units == 'kg':
        return fuel_heat_value*original_value*boiler_efficiency

def compute_cost_per_kJ(fuel, cost, cost_units, boiler_type):
    if cost_units == '$/kWh':
        return cost/3600
    if fuel == 'Electricidad':
        raise ValueError('El costo de la electricidad solo puede expresarse en $/kWh')
    if fuel == 'Propano':
        fuel_density = 1.83
        if boiler_type == 'Condensación líquido' or boiler_type == 'Condensación vapor':
            fuel_heat_value = 50426
        else:
            fuel_heat_value = 46367
    if fuel == 'Butano':
        fuel_density = 2.52
        if boiler_type == 'Condensación líquido' or boiler_type == 'Condensación vapor':
            fuel_heat_value = 49652
        else:
            fuel_heat_value = 45765
    if fuel == 'Biomasa':
        fuel_density = 800
        if boiler_type == 'Condensación líquido' or boiler_type == 'Condensación vapor':
            fuel_heat_value = 19000
        else:
            fuel_heat_value = 17000
    if fuel == 'Diesel':
        fuel_density = 850
        if boiler_type == 'Condensación líquido' or boiler_type == 'Condensación vapor':
            fuel_heat_value = 45600
        else:
            fuel_heat_value = 42600
    if cost_units == '$/L':
        kJ_per_litre = fuel_density*fuel_heat_value/1000
        litres_per_kJ = 1/kJ_per_litre
        return cost*litres_per_kJ
    if cost_units == '$/m3':
        kJ_per_m3 = fuel_density*fuel_heat_value
        m3_per_kJ = 1/kJ_per_m3
        return cost*m3_per_kJ
    if cost_units == '$/kg':
        kg_per_kJ = 1/fuel_heat_value
        return cost*kg_per_kJ

def time_from_string(string):
    List = string.split(':')
    time = int(List[0])
    if List[1] == '30':
        time = np.round(time + 0.5, 1)
    return time

def week_day_from_day(day):
    week_day = day%7
    if week_day == 1:
        return 'Monday'
    if week_day == 2:
        return 'Tuesday'
    if week_day == 3:
        return 'Wednesday'
    if week_day == 4:
        return 'Thursday'
    if week_day == 5:
        return 'Friday'
    if week_day == 6:
        return 'Saturday'
    if week_day == 0:
        return 'Sunday'

def month_from_day(day):
    if day <= 31:
        return 'January'
    if day >= 32 and day <= 59:
        return 'February'
    if day >= 60 and day <= 90:
        return 'March'
    if day >= 91 and day <= 120:
        return 'April'
    if day >= 121 and day <= 151:
        return 'May'
    if day >= 152 and day <= 181:
        return 'June'
    if day >= 182 and day <= 212:
        return 'July'
    if day >= 213 and day <= 243:
        return 'August'
    if day >= 244 and day <= 273:
        return 'September'
    if day >= 274 and day <= 304:
        return 'October'
    if day >= 305 and day <= 334:
        return 'November'
    if day >= 335:
        return 'December'
    
def compute_VAN(initial_investment, yearly_savings, discount_factor, inflation, tax_rate, initial_dep_value, dep_time):
    if discount_factor > 1:
        discount_factor = discount_factor/100
    if inflation > 1:
        inflation = inflation/100
    if tax_rate > 1:
        tax_rate = tax_rate/100
    VAN = -initial_investment
    for year in range(1,21):
        if year <= dep_time:
            utility = yearly_savings - initial_dep_value/dep_time
        else:
            utility = yearly_savings
        utility = (1 - tax_rate)*utility
        if year <= dep_time:
            utility = utility + initial_dep_value/dep_time
        VAN = VAN + utility/( (1 + discount_factor)**year )
    return VAN

def VAN_squared(discount_factor, initial_investment, yearly_savings, inflation, tax_rate, initial_dep_value, dep_time):
    return (compute_VAN(initial_investment, yearly_savings, discount_factor, inflation, tax_rate, initial_dep_value, dep_time))**2
    
def compute_TIR(initial_investment, yearly_savings, inflation, tax_rate, initial_dep_value, dep_time):
    discount_factor_bounds = Bounds(lb = 0, ub = np.inf)
    solution = minimize(VAN_squared, x0 = 7, args = (initial_investment, yearly_savings, inflation, tax_rate, initial_dep_value, dep_time), bounds = discount_factor_bounds)
    return solution.x[0]


