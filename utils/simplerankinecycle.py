from CoolProp.CoolProp import PropsSI
from CoolProp.Plots.SimpleCycles import StateContainer

class RankineCycle:
    """
    Class to compute and store the thermodynamic states of a Rankine cycle.
    This only works with saturated phase only
    """

    def __init__(self, fluid, eta_exp, eta_pum):
        """
        Initialize the RankineCycle class.

        Parameters:
        -----------
        fluid : str
            Working fluid (e.g., 'Water').
        eta_exp : float
            Isentropic efficiency of the expander (turbine).
        eta_pum : float
            Isentropic efficiency of the pump.
        """
        self.fluid = fluid
        self.eta_exp = eta_exp
        self.eta_pum = eta_pum
        self.cycle_states = StateContainer(unit_system="SI")

    def compute_states(self, p2, q2, p0, q0):
        """
        Compute the thermodynamic states of the Rankine cycle.

        Parameters:
        -----------
        p2 : float
            Turbine inlet pressure (Pa).
        q2 : float
            Turbine inlet quality (1.0 for saturated vapor).
        p0 : float
            Condenser outlet pressure (Pa).
        q0 : float
            Condenser outlet quality (0.0 for saturated liquid).
        """
        # State 2: Saturated vapor entering the turbine
        h2 = PropsSI("H", "P", p2, "Q", q2, self.fluid)
        s2 = PropsSI("S", "P", p2, "Q", q2, self.fluid)
        T2 = PropsSI("T", "P", p2, "Q", q2, self.fluid)
        self._store_state(2, p2, T2, h2, s2)

        # State 3: Expanded gas exiting the turbine
        p3 = p0
        h3_isentropic = PropsSI("H", "P", p3, "S", s2, self.fluid)
        h3 = h2 - self.eta_exp * (h2 - h3_isentropic)
        T3 = PropsSI("T", "P", p3, "H", h3, self.fluid)
        s3 = PropsSI("S", "P", p3, "H", h3, self.fluid)
        self._store_state(3, p3, T3, h3, s3)

        # State 0: Saturated liquid exiting the condenser
        h0 = PropsSI("H", "P", p0, "Q", q0, self.fluid)
        s0 = PropsSI("S", "P", p0, "Q", q0, self.fluid)
        T0 = PropsSI("T", "P", p0, "Q", q0, self.fluid)
        self._store_state(0, p0, T0, h0, s0)

        # State 1: Compressed liquid after the pump
        p1 = p2
        h1_isentropic = PropsSI("H", "P", p1, "S", s0, self.fluid)
        h1 = h0 + (h1_isentropic - h0) / self.eta_pum
        T1 = PropsSI("T", "P", p1, "H", h1, self.fluid)
        s1 = PropsSI("S", "P", p1, "H", h1, self.fluid)
        self._store_state(1, p1, T1, h1, s1)

    def _store_state(self, state_id, p, T, h, s):
        """
        Store the thermodynamic state in the cycle_states container.

        Parameters:
        -----------
        state_id : int
            State identifier (e.g., 0, 1, 2, 3).
        p : float
            Pressure (Pa).
        T : float
            Temperature (K).
        h : float
            Enthalpy (J/kg).
        s : float
            Entropy (J/(kg·K)).
        """
        self.cycle_states[state_id, 'P'] = p
        self.cycle_states[state_id, 'T'] = T
        self.cycle_states[state_id, 'H'] = h
        self.cycle_states[state_id, 'S'] = s

    def get_states(self):
        """
        Retrieve the computed cycle states.

        Returns:
        --------
        StateContainer
            Container with all computed thermodynamic states.
        """
        return self.cycle_states
