import matplotlib.pyplot as plt

def heat_exchanger(stream_hot, stream_cold, T0):
    """
    Calculate exergy destruction in a heat exchanger.
    Assumptions:
        - Steady-state operation
        - No heat losses to the environment
        - Negligible kinetic and potential energy differences

    Parameters:
        stream_hot (dict): {
            "m_dot": float,      # mass flow rate (kg/s)
            "h": {"in": float, "out": float},   # enthalpy (J/kg)
            "s": {"in": float, "out": float}    # entropy (J/kg.K)
        }
        stream_cold (dict): {
            "m_dot": float,
            "h": {"in": float, "out": float},
            "s": {"in": float, "out": float}
        }
        T0 (float): Ambient (dead state) temperature (K)

    Returns:
        float: Exergy destruction (W)
    """
    # Flow exergy change for hot stream
    ex_hot = stream_hot["m_dot"] * (
        (stream_hot["h"]["in"] - stream_hot["h"]["out"])
        - T0 * (stream_hot["s"]["in"] - stream_hot["s"]["out"])
    )

    # Flow exergy change for cold stream
    ex_cold = stream_cold["m_dot"] * (
        (stream_cold["h"]["in"] - stream_cold["h"]["out"])
        - T0 * (stream_cold["s"]["in"] - stream_cold["s"]["out"])
    )

    # Exergy destruction
    return ex_hot + ex_cold

def work_devices(device_stream, specific_work, gravity, T0):
    """
    Calculate specific exergy destruction for work-producing or work-consuming devices (e.g., turbines, compressors, pumps).
    Assumptions:
        - Steady-state operation
        - No heat losses to the environment
        - Single inlet, single outlet
    Parameters:
        device_stream (dict): {
            "h": {"in": float, "out": float},      # enthalpy (J/kg)
            "s": {"in": float, "out": float},      # entropy (J/kg.K)
            "vel": {"in": float, "out": float},    # velocity (m/s)
            "pos": {"in": float, "out": float}     # position (m)
        }
        specific_work (float): Specific work input/output (J/kg)
        gravity (float): Gravitational acceleration (m/s^2)
        T0 (float): Ambient (dead state) temperature (K)

    Returns:
        float: Exergy destruction (W/kg)
    """
    h_in = device_stream["h"]["in"]
    h_out = device_stream["h"]["out"]
    s_in = device_stream["s"]["in"]
    s_out = device_stream["s"]["out"]
    vel_in = device_stream["vel"]["in"]
    vel_out = device_stream["vel"]["out"]
    pos_in = device_stream["pos"]["in"]
    pos_out = device_stream["pos"]["out"]

    ex_flow = (h_in - h_out) - T0 * (s_in - s_out) + (vel_in**2 - vel_out**2) / 2 + gravity * (pos_in - pos_out)
    return -specific_work + ex_flow

def cycle_exergy_efficiency(net_work_output, heat_source_stream, T0):
    """
    Calculate the overall cycle exergy efficiency.
    Assumptions:
        - Steady-state operation
        - No heat losses to the environment
        - Single heat source stream
        - Negligible kinetic and potential energy differences

    Parameters:
        net_work_output (float): Net work output of the cycle (W)
        heat_source_stream (dict): {
            "h": {"in": float, "out": float},   # enthalpy (J/kg)
            "s": {"in": float, "out": float}    # entropy (J/kg.K)
        }
        T0 (float): Ambient (dead state) temperature (K)

    Returns:
        float: Exergy efficiency of the cycle (dimensionless, between 0 and 1)
    """
    h_in = heat_source_stream["h"]["in"]
    h_out = heat_source_stream["h"]["out"]
    s_in = heat_source_stream["s"]["in"]
    s_out = heat_source_stream["s"]["out"]
    ex_flow = h_in - h_out - T0 * (s_in - s_out)
    return net_work_output / ex_flow

def plot_exergy_bar_chart(exergy_dict, title="Exergy Destruction by Component"):
    """
    Plots a bar chart of exergy destruction for each component.

    Parameters:
        exergy_dict (dict): Keys are component names, values are exergy values (in kW).
        title (str): Title of the chart.
    """
    components = list(exergy_dict.keys())
    exergy_values = [val / 1000 for val in exergy_dict.values()]

    plt.figure(figsize=(10, 6))
    bars = plt.bar(components, exergy_values, color='skyblue')
    plt.xlabel("Component")
    plt.ylabel("Exergy Destruction (kW)")
    plt.title(title)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    # Annotate bars with values
    for bar, value in zip(bars, exergy_values):
        plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{value:.2f}", 
                    ha='center', va='bottom', fontsize=9)

    plt.show()

def plot_exergy_pie_chart(exergy_dict, title="Exergy Destruction Distribution"):
    """
    Plots a pie chart of exergy destruction distribution by component.

    Parameters:
        exergy_dict (dict): Keys are component names, values are exergy values (in kW).
        title (str): Title of the chart.
    """
    components = list(exergy_dict.keys())
    exergy_values = [val / 1000 for val in exergy_dict.values()]

    plt.figure(figsize=(8, 8))
    plt.pie(
        exergy_values,
        labels=components,
        autopct='%1.1f%%',
        startangle=140,
        colors=plt.cm.Paired.colors
    )
    plt.title(title)
    plt.tight_layout()
    plt.show()

# # # Example usage with dummy data
# exergy_data = {
#     "Boiler": 120.5,
#     "Condenser": 85.3,
#     "Turbine": 210.7,
#     "Pump": 45.2,
#     "Heat Exchanger": 60.9
# }

# plot_exergy_bar_chart(exergy_data, title="ORC Exergy Destruction by Component")

# plot_exergy_pie_chart(exergy_data, title="ORC Exergy Destruction Distribution")
