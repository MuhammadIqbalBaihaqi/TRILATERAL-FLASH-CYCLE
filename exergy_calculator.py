

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

