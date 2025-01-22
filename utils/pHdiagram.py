# # p-H diagram of working fluids
# # Contributor: Matahari Arsyabil Muhammad Choiri, Sherryn, Muhammad Jati, Muhammad Zelot Zoha
# import math
# from numpy import *
# from matplotlib.pyplot import *

# import CoolProp
# from CoolProp.CoolProp import PropsSI

# def pHdiagram(Fluid):
#     """
#     Plot the saturation curve of a given fluid in p-H coordinates.

#     Parameters:
#     Fluid (str): Name of the fluid (e.g., 'propane').

#     """
#     # Get critical and triple point temperatures
#     T_crit = PropsSI('Tcrit', Fluid)
#     T_triple = PropsSI('Ttriple', Fluid)

#     # Generate temperature range
#     T = linspace(T_triple, T_crit, 1000)

#     # Initialize lists for pressure and enthalpy values
#     p_0 = []  # Pressure for saturated liquid (Q=0)
#     h_0 = []  # Enthalpy for saturated liquid (Q=0)
#     p_1 = []  # Pressure for saturated vapor (Q=1)
#     h_1 = []  # Enthalpy for saturated vapor (Q=1)

#     # Calculate properties at each temperature in the range
#     for temp in T:
#         p_0.append(PropsSI('P', 'T', temp, 'Q', 0, Fluid))
#         h_0.append(PropsSI('H', 'T', temp, 'Q', 0, Fluid))
#         p_1.append(PropsSI('P', 'T', temp, 'Q', 1, Fluid))
#         h_1.append(PropsSI('H', 'T', temp, 'Q', 1, Fluid))

#     # Convert pressure to MPa and enthalpy to kJ/kg for better readability
#     p_0 = [p / 1e6 for p in p_0]
#     p_1 = [p / 1e6 for p in p_1]
#     h_0 = [h / 1000 for h in h_0]
#     h_1 = [h / 1000 for h in h_1]

#     # Combine data for a closed-loop plot
#     h_combined = h_0 + h_1[::-1]
#     p_combined = concatenate([p_0, p_1[::-1]])

#     # Plot the saturation graph in p-H coordinates with logarithmic pressure scale
#     figure()
#     plot(h_combined, p_combined, 'k')
#     xlabel('$h$ (kJ/kg)', fontsize=12)
#     ylabel('$p$ (MPa)', fontsize=12)

#     # Set X and Y limits dynamically based on calculated values
#     x_min, x_max = min(h_combined), max(h_combined)
#     x_min = math.floor(x_min/1000)*1000
#     x_max = math.ceil(x_max/1000)*1000

#     y_min_real, y_max_real = min(p_combined), max(p_combined)
#     y_min_power=math.floor(math.log10(y_min_real))
#     y_min=10**y_min_power
#     y_max_power=math.ceil(math.log10(y_max_real))
#     y_max=10**y_max_power

#     if Fluid.upper()=='MD4M':
#       y_max=10
#     elif Fluid.lower()=='toluene':
#       y_max=100
#     xlim(x_min, x_max)  # Sedikit memperbesar batas sumbu X
#     ylim(y_min, y_max)  # Sedikit memperbesar batas sumbu Y
#     yscale('log')  # Set the pressure axis to logarithmic scale

#     title(Fluid)
#     grid(False)

#     # Add legend and adjust layout
#     tight_layout()
#     show()


import math
from numpy import *
from matplotlib.pyplot import *

import CoolProp
from CoolProp.CoolProp import PropsSI

def pHdiagram(Fluid, extracted_ph=None):
    """
    Plot the saturation curve of a given fluid in p-H coordinates and optionally overlay extracted p-H data.

    Parameters:
    Fluid (str): Name of the fluid (e.g., 'propane').
    extracted_ph (array-like): Optional p-H data (array of [p, h] values) to overlay on the diagram.

    """
    # Get critical and triple point temperatures
    T_crit = PropsSI('Tcrit', Fluid)
    T_triple = PropsSI('Ttriple', Fluid)

    # Generate temperature range
    T = linspace(T_triple, T_crit, 1000)

    # Initialize lists for pressure and enthalpy values
    p_0 = []  # Pressure for saturated liquid (Q=0)
    h_0 = []  # Enthalpy for saturated liquid (Q=0)
    p_1 = []  # Pressure for saturated vapor (Q=1)
    h_1 = []  # Enthalpy for saturated vapor (Q=1)

    # Calculate properties at each temperature in the range
    for temp in T:
        p_0.append(PropsSI('P', 'T', temp, 'Q', 0, Fluid))
        h_0.append(PropsSI('H', 'T', temp, 'Q', 0, Fluid))
        p_1.append(PropsSI('P', 'T', temp, 'Q', 1, Fluid))
        h_1.append(PropsSI('H', 'T', temp, 'Q', 1, Fluid))

    # Convert pressure to MPa and enthalpy to kJ/kg for better readability
    p_0 = [p / 1e6 for p in p_0]
    p_1 = [p / 1e6 for p in p_1]
    h_0 = [h / 1000 for h in h_0]
    h_1 = [h / 1000 for h in h_1]

    # Combine data for a closed-loop plot
    h_combined = h_0 + h_1[::-1]
    p_combined = concatenate([p_0, p_1[::-1]])

    # Plot the saturation graph in p-H coordinates with logarithmic pressure scale
    figure()
    plot(h_combined, p_combined, 'k', label='Saturation Curve')

    # If extracted p-H data is provided, overlay it on the plot
    if extracted_ph is not None:
        p_extracted, h_extracted = zip(*extracted_ph)
        p_extracted = [p for p in p_extracted]  
        h_extracted = [h for h in h_extracted]  
        plot(h_extracted, p_extracted, color='blue', linestyle='solid', linewidth=1, label='Process Cycle')

    xlabel('$h$ (kJ/kg)', fontsize=12)
    ylabel('$p$ (MPa)', fontsize=12)

    # Set X and Y limits dynamically based on calculated values
    x_min, x_max = min(h_combined), max(h_combined)
    x_min = math.floor(x_min/1000)*1000
    x_max = math.ceil(x_max/1000)*1000

    y_min_real, y_max_real = min(p_combined), max(p_combined)
    y_min_power = math.floor(math.log10(y_min_real))
    y_min = 10**y_min_power
    y_max_power = math.ceil(math.log10(y_max_real))
    y_max = 10**y_max_power

    if Fluid.upper() == 'MD4M':
        y_max = 10
    elif Fluid.lower() == 'toluene':
        y_max = 100

    xlim(x_min, x_max)  # Sedikit memperbesar batas sumbu X
    ylim(y_min, y_max)  # Sedikit memperbesar batas sumbu Y
    yscale('log')  # Set the pressure axis to logarithmic scale

    title(Fluid)
    grid(False)

    # Add legend and adjust layout
    legend()
    tight_layout()
    show()

# Example usage:
# pHdiagram('Water', extracted_ph=[(101325, 500000), (200000, 700000), (300000, 900000)])
