# # T-s diagram of working fluids
# # Contributor: Matahari Arsyabil Muhammad Choiri, Sherryn, Muhammad Jati, Muhammad Zelot Zoha
# from numpy import *
# from matplotlib.pyplot import *
# import math
# import CoolProp
# from CoolProp.CoolProp import PropsSI

# def TSdiagram(Fluid):
#     """
#     Plot the saturation curve of a given fluid in T-s coordinates.

#     Parameters:
#     Fluid (str): Name of the fluid (e.g., 'propane').

#     """
#     # Get critical and triple point temperatures
#     T_crit = PropsSI('Tcrit', Fluid)
#     T_triple = PropsSI('Ttriple', Fluid)

#     # Generate temperature range
#     T = linspace(T_triple, T_crit, 1000)

#     # Initialize lists for entropy values
#     s_1 = []  # Entropy for saturated vapor (Q=1)
#     s_0 = []  # Entropy for saturated liquid (Q=0)

#     # Calculate properties at each temperature in the range
#     for temp in T:
#         s_1.append(PropsSI('S', 'T', temp, 'Q', 1, Fluid))
#         s_0.append(PropsSI('S', 'T', temp, 'Q', 0, Fluid))

#     # Convert entropy from J/(kg⋅K) to kJ/(kg⋅K) for better readability
#     s_0 = [s / 1000 for s in s_0]
#     s_1 = [s / 1000 for s in s_1]

#     # Combined data for a closed-loop plot
#     s_combined = s_0 + s_1[::-1]
#     T_combined = concatenate([T, T[::-1]])

#     # Plot the saturation graph in T-s coordinates
#     figure()
#     plot(s_combined, T_combined, 'k')
#     xlabel('$s$ (kJ/(kg⋅K)', fontsize=12)
#     ylabel('$T$ (K)', fontsize=12)
#     x_min, x_max = min(s_combined), max(s_combined)
#     y_min, y_max = min(T_combined), max(T_combined)

#     x_min = math.floor(x_min)
#     x_max = math.ceil(x_max)

#     y_min = math.floor(y_min/100)*100
#     y_max = math.ceil(y_max/100)*100

#     xlim(x_min, x_max)  # Sedikit memperbesar batas sumbu X
#     ylim(y_min, y_max)  # Sedikit memperbesar batas sumbu Y
#     title(Fluid)
#     grid(False)

#     # Adjust layout to prevent axis labels from being clipped
#     tight_layout()


from numpy import *
from matplotlib.pyplot import *
import math
import CoolProp
from CoolProp.CoolProp import PropsSI

def TSdiagram(Fluid, extracted_ts=None):
    """
    Plot the saturation curve of a given fluid in T-s coordinates.

    Parameters:
    Fluid (str): Name of the fluid (e.g., 'propane').
    extracted_ts (array-like): Optional T-s data (array of [T, s] values) to overlay on the diagram.

    """
    # Get critical and triple point temperatures
    T_crit = PropsSI('Tcrit', Fluid)
    T_triple = PropsSI('Ttriple', Fluid)

    # Generate temperature range
    T = linspace(T_triple, T_crit, 1000)

    # Initialize lists for entropy values
    s_1 = []  # Entropy for saturated vapor (Q=1)
    s_0 = []  # Entropy for saturated liquid (Q=0)

    # Calculate properties at each temperature in the range
    for temp in T:
        s_1.append(PropsSI('S', 'T', temp, 'Q', 1, Fluid) / 1000)  # Convert to kJ/kg.K
        s_0.append(PropsSI('S', 'T', temp, 'Q', 0, Fluid) / 1000)  # Convert to kJ/kg.K

    # Combined data for a closed-loop plot
    s_combined = s_0 + s_1[::-1]
    T_combined = concatenate([T, T[::-1]])

    # Plot the saturation graph in T-s coordinates
    figure()
    plot(s_combined, T_combined, 'k', label='Saturation Curve')

    # If extracted T-s data is provided, overlay it on the plot
    if extracted_ts is not None:
        T_extracted, s_extracted = zip(*extracted_ts)
        s_extracted = [s for s in s_extracted]  
        plot(s_extracted, T_extracted, color='blue', linestyle='solid', linewidth=1, label='Process Cycle')


    xlabel('$s$ (kJ/(kg⋅K))', fontsize=12)
    ylabel('$T$ (K)', fontsize=12)
    x_min, x_max = min(s_combined), max(s_combined)
    y_min, y_max = min(T_combined), max(T_combined)

    x_min = math.floor(x_min)
    x_max = math.ceil(x_max)

    y_min = math.floor(y_min/100)*100
    y_max = math.ceil(y_max/100)*100

    xlim(x_min, x_max)  # Sedikit memperbesar batas sumbu X
    ylim(y_min, y_max)  # Sedikit memperbesar batas sumbu Y
    title(Fluid)
    grid(False)
    legend()

    # Adjust layout to prevent axis labels from being clipped
    tight_layout()

# Example usage:
# TSdiagram('Water', extracted_ts=[(300, 500), (350, 600), (400, 700)])
