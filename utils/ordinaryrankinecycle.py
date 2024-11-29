from CoolProp.CoolProp import PropsSI
from CoolProp.Plots.SimpleCycles import StateContainer


class ReheatRankineCycle:
    """
    Class to compute and store thermodynamic states for a Rankine cycle with reheating.
    """
    STATECOUNT = 6
    def __init__(self, fluid, eta_turbine_high, eta_turbine_low, eta_pump):
        """
        Initialize the Rankine cycle parameters.

        Parameters:
        -----------
        fluid : str
            Working fluid (e.g., 'Water').
        eta_turbine_high : float
            Isentropic efficiency of the high-pressure turbine.
        eta_turbine_low : float
            Isentropic efficiency of the low-pressure turbine.
        eta_pump : float
            Isentropic efficiency of the pump.
        """
        self.fluid = fluid
        self.eta_turbine_high = eta_turbine_high
        self.eta_turbine_low = eta_turbine_low
        self.eta_pump = eta_pump
        self.cycle_states = StateContainer(unit_system="SI")

    def compute_cycle(self, p_high_in, T_high_in, p_high_out, T_reheat_in, p_cond, q_cond):
        """
        Compute the thermodynamic states for the Rankine cycle with reheating.

        Parameters:
        -----------
        p_high_in : float
            Pressure before the high-pressure turbine (Pa).
        T_high_in : float
            Temperature before the high-pressure turbine (K).
        p_high_out : float
            Pressure after the high-pressure turbine (Pa).
        T_reheat_in : float
            Temperature after reheating and before the low-pressure turbine (K).
        p_cond : float
            Condenser (isobar) pressure (Pa).
        q_cond : float
            Quality after the condenser (0.0 for saturated liquid).
        """
        # State 1: Saturated liquid at condenser pressure
        self._store_state(1, p_cond, PropsSI("T", "P", p_cond, "Q", q_cond, self.fluid),
                          PropsSI("H", "P", p_cond, "Q", q_cond, self.fluid),
                          PropsSI("S", "P", p_cond, "Q", q_cond, self.fluid))

        # State 2: Compressed liquid after the pump
        p2 = p_high_in
        h2_isentropic = PropsSI("H", "P", p2, "S", self.cycle_states[1, 'S'], self.fluid)
        h2 = self.cycle_states[1, 'H'] + (h2_isentropic - self.cycle_states[1, 'H']) / self.eta_pump
        self._store_state(2, p2, PropsSI("T", "P", p2, "H", h2, self.fluid), h2,
                          PropsSI("S", "P", p2, "H", h2, self.fluid))

        # State 3: High-pressure turbine inlet
        self._store_state(3, p_high_in, T_high_in,
                          PropsSI("H", "P", p_high_in, "T", T_high_in, self.fluid),
                          PropsSI("S", "P", p_high_in, "T", T_high_in, self.fluid))

        # State 4: After high-pressure turbine
        h4_isentropic = PropsSI("H", "P", p_high_out, "S", self.cycle_states[3, 'S'], self.fluid)
        h4 = self.cycle_states[3, 'H'] - self.eta_turbine_high * (self.cycle_states[3, 'H'] - h4_isentropic)
        self._store_state(4, p_high_out, PropsSI("T", "P", p_high_out, "H", h4, self.fluid), h4,
                          PropsSI("S", "P", p_high_out, "H", h4, self.fluid))

        # State 5: Low-pressure turbine inlet (reheat stage)
        self._store_state(5, p_high_out, T_reheat_in,
                          PropsSI("H", "P", p_high_out, "T", T_reheat_in, self.fluid),
                          PropsSI("S", "P", p_high_out, "T", T_reheat_in, self.fluid))

        # State 6: After low-pressure turbine
        h6_isentropic = PropsSI("H", "P", p_cond, "S", self.cycle_states[5, 'S'], self.fluid)
        h6 = self.cycle_states[5, 'H'] - self.eta_turbine_low * (self.cycle_states[5, 'H'] - h6_isentropic)
        self._store_state(6, p_cond, PropsSI("T", "P", p_cond, "H", h6, self.fluid), h6,
                          PropsSI("S", "P", p_cond, "H", h6, self.fluid))

    def _store_state(self, state_id, p, T, h, s):
        """
        Store a thermodynamic state in the state container.

        Parameters:
        -----------
        state_id : int
            State identifier.
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
