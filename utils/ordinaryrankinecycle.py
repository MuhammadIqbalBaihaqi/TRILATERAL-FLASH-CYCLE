from CoolProp.CoolProp import PropsSI
from CoolProp.Plots.SimpleCycles import StateContainer


from CoolProp.CoolProp import PropsSI
from CoolProp.Plots.SimpleCycles import StateContainer, BaseCycle

class ReheatRankineCycle(BaseCycle):
    """
    Class to compute and store thermodynamic states for a Rankine cycle with reheating.
    """
    STATECOUNT = 6
    STATECHANGE = [
        lambda inp: BaseCycle.state_change(inp, 'S', 'P', 0, ty1='log', ty2='log'),  # Pumping process
        lambda inp: BaseCycle.state_change(inp, 'H', 'P', 1, ty1='lin', ty2='lin'),  # Heat addition (boiler)
        lambda inp: BaseCycle.state_change(inp, 'H', 'P', 2, ty1='log', ty2='log'),  # High-pressure turbine expansion
        lambda inp: BaseCycle.state_change(inp, 'H', 'P', 3, ty1='lin', ty2='lin'),  # Reheating
        lambda inp: BaseCycle.state_change(inp, 'H', 'P', 4, ty1='log', ty2='log'),  # Low-pressure turbine expansion
        lambda inp: BaseCycle.state_change(inp, 'H', 'P', 5, ty1='lin', ty2='lin')   # Condenser heat removal
    ]

    def __init__(self, fluid_ref, eta_turbine_high, eta_turbine_low, eta_pump, graph_type="TS", **kwargs):
        """
        Initialize the Rankine cycle with reheating.

        Parameters:
        -----------
        fluid_ref : str
            Working fluid (e.g., 'Water').
        eta_turbine_high : float
            Isentropic efficiency of the high-pressure turbine.
        eta_turbine_low : float
            Isentropic efficiency of the low-pressure turbine.
        eta_pump : float
            Isentropic efficiency of the pump.
        """
        self.eta_turbine_high = eta_turbine_high
        self.eta_turbine_low = eta_turbine_low
        self.eta_pump = eta_pump
        self.fluid_ref = fluid  # Align with the inherited structure
        super().__init__(fluid, graph_type, **kwargs)
        # self.eta_turbine_high = eta_turbine_high
        # self.eta_turbine_low = eta_turbine_low
        # self.eta_pump = eta_pump
        # super().__init__(fluid_ref, graph_type, **kwargs)

    def simple_solve(self, p_high_in, T_high_in, p_high_out, T_reheat_in, p_cond, q_cond, SI=True):
        """
        Solve the thermodynamic states for a Rankine cycle with reheating.

        Parameters:
        -----------
        p_high_in : float
            High-pressure turbine inlet pressure (Pa).
        T_high_in : float
            High-pressure turbine inlet temperature (K).
        p_high_out : float
            High-pressure turbine outlet pressure (Pa).
        T_reheat_in : float
            Low-pressure turbine inlet temperature (K).
        p_cond : float
            Condenser pressure (Pa).
        q_cond : float
            Quality after the condenser (e.g., 0.0 for saturated liquid).
        """
        cycle_states = StateContainer(unit_system=self._system)

        # State 0: Saturated liquid at condenser pressure
        T0 = PropsSI("T", "P", p_cond, "Q", q_cond, self.fluid_ref)
        h0 = PropsSI("H", "P", p_cond, "Q", q_cond, self.fluid_ref)
        s0 = PropsSI("S", "P", p_cond, "Q", q_cond, self.fluid_ref)
        cycle_states[0, 'P'], cycle_states[0, 'T'], cycle_states[0, 'H'], cycle_states[0, 'S'] = p_cond, T0, h0, s0

        # State 1: Pump outlet
        h1_isentropic = PropsSI("H", "P", p_high_in, "S", s0, self.fluid_ref)
        h1 = h0 + (h1_isentropic - h0) / self.eta_pump
        self.state.update(CoolProp.HmassP_INPUTS, h1, p_high_in)
        T1 = self.state.T()
        s1 = self.state.smass()
        cycle_states[1, 'P'], cycle_states[1, 'T'], cycle_states[1, 'H'], cycle_states[1, 'S'] = p_high_in, T1, h1, s1

        # State 2: High-pressure turbine inlet
        self.state.update(CoolProp.PT_INPUTS, p_high_in, T_high_in)
        h2 = self.state.hmass()
        s2 = self.state.smass()
        cycle_states[2, 'P'], cycle_states[2, 'T'], cycle_states[2, 'H'], cycle_states[2, 'S'] = p_high_in, T_high_in, h2, s2

        # State 3: High-pressure turbine outlet
        h3_isentropic = PropsSI("H", "P", p_high_out, "S", s2, self.fluid_ref)
        h3 = h2 - self.eta_turbine_high * (h2 - h3_isentropic)
        self.state.update(CoolProp.HmassP_INPUTS, h3, p_high_out)
        T3 = self.state.T()
        s3 = self.state.smass()
        cycle_states[3, 'P'], cycle_states[3, 'T'], cycle_states[3, 'H'], cycle_states[3, 'S'] = p_high_out, T3, h3, s3

        # State 4: Reheat stage
        self.state.update(CoolProp.PT_INPUTS, p_high_out, T_reheat_in)
        h4 = self.state.hmass()
        s4 = self.state.smass()
        cycle_states[4, 'P'], cycle_states[4, 'T'], cycle_states[4, 'H'], cycle_states[4, 'S'] = p_high_out, T_reheat_in, h4, s4

        # State 5: Low-pressure turbine outlet
        h5_isentropic = PropsSI("H", "P", p_cond, "S", s4, self.fluid_ref)
        h5 = h4 - self.eta_turbine_low * (h4 - h5_isentropic)
        self.state.update(CoolProp.HmassP_INPUTS, h5, p_cond)
        T5 = self.state.T()
        s5 = self.state.smass()
        cycle_states[5, 'P'], cycle_states[5, 'T'], cycle_states[5, 'H'], cycle_states[5, 'S'] = p_cond, T5, h5, s5

        self.cycle_states = cycle_states
        self.fill_states()

    def eta_thermal(self):
        """
        Calculate the thermal efficiency of the cycle.

        Returns:
        --------
        float
            Thermal efficiency.
        """
        w_net = self.cycle_states[2].H - self.cycle_states[3].H + self.cycle_states[4].H - self.cycle_states[5].H
        q_boiler = self.cycle_states[2].H - self.cycle_states[1].H + self.cycle_states[4].H - self.cycle_states[3].H
        return w_net / q_boiler



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
