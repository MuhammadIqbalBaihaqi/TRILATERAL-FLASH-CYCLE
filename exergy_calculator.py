

def heat_exchanger(stream_hot, stream_cold, T0):
    """
    Calculate exergy destruction in a heat exchanger.

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

