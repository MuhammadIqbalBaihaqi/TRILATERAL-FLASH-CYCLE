from CoolProp.CoolProp import PropsSI
import CoolProp
from CoolProp.Plots.SimpleCycles import StateContainer, BaseCycle
import numpy as np

class ORCSuperheat(BaseCycle):
    """
    Class to compute and store thermodynamic states for a Rankine cycle with superheating, a desuperheater, and a pump.
    """
    STATECOUNT = 7  # 7 states for the Rankine cycle based on the corrected instructions
    STATECHANGE = [
        lambda inp: BaseCycle.state_change(inp, 'H', 'P', 0, ty1='lin', ty2='lin'),  # Isobaric heating (Heater)
        lambda inp: BaseCycle.state_change(inp, 'H', 'P', 1, ty1='lin', ty2='lin'),  # Isobaric heating (Superheater)
        lambda inp: BaseCycle.state_change(inp, 'H', 'P', 2, ty1='log', ty2='log'),  # Isentropic expansion (Turbine)
        lambda inp: BaseCycle.state_change(inp, 'H', 'P', 3, ty1='lin', ty2='lin'),  # Isobaric cooling (Desuperheater) 
        lambda inp: BaseCycle.state_change(inp, 'H', 'P', 4, ty1='lin', ty2='lin'),  # Isobaric cooling (Condenser)
        lambda inp: BaseCycle.state_change(inp, 'S', 'P', 5, ty1='log', ty2='log'),  # Isentropic compression (Pump)
        lambda inp: BaseCycle.state_change(inp, 'H', 'P', 6, ty1='lin', ty2='lin'),  # Isobaric heating (Preheater)
    ]

    def __init__(self, fluid, eta_turbine, eta_pump, m_dot_hs, graph_type="TS", **kwargs):
        """
        Initialize the Rankine cycle with superheating, desuperheating, and a pump.

        Parameters:
        -----------
        fluid : str
            Working fluid (e.g., 'Water').
        eta_turbine : float
            Isentropic efficiency of the turbine.
        eta_pump : float
            Isentropic efficiency of the pump.
        """
        self.eta_turbine = eta_turbine
        self.eta_pump = eta_pump
        self.fluid_ref = fluid  # Align with the inherited structure
        self.m_dot_hs = m_dot_hs
        super().__init__(fluid, graph_type, **kwargs)

    def simple_solve(self, T_s_in, pp_1sup, T_s_mid1, pp_0, T_cs_mid1, pp_cs_mid1, T_cs_in, SI=True):
        """
        Solve the thermodynamic states for the Rankine cycle with superheating, desuperheating, and a pump.

        Parameters:
        -----------
        T_s_in : float
            Heat source temperature (K).
        pp_1sup : float
            Superheater pinch point (K).
        T_s_mid1 : float
            Middle source temperature (K).
        pp_0 : float
            Middle source pinch point (K).
        T_cs_mid1 : float
            Middle cooling source temperature (K).
        T_cs_in : float
            Inlet cooling source temperature (K).
        pp_cs_mid1 : float
            Middle cooling source pinch point (K).
        """
        cycle_states = StateContainer(unit_system=self._system)

        # State 0: After preheater (saturated liquid, high pressure)
        T_0 = T_s_mid1 - pp_0
        p_0 = PropsSI('P', 'Q', 0.0, 'T', T_0,self.fluid_ref)
        h_0 = PropsSI('H', 'Q', 0.0, 'T', T_0, self.fluid_ref)
        s_0 = PropsSI('S', 'Q', 0.0, 'T', T_0, self.fluid_ref)
        cycle_states[0, 'P'], cycle_states[0, 'T'], cycle_states[0, 'H'], cycle_states[0, 'S'] = p_0, T_0, h_0, s_0

        # State 1: After heater (saturated vapor, high pressure)
        T_1 = T_0
        p_1 = PropsSI('P', 'Q', 1.0, 'T', T_1, self.fluid_ref)
        h_1 = PropsSI('H', 'Q', 1.0, 'T', T_1, self.fluid_ref)
        s_1 = PropsSI('S', 'Q', 1.0, 'T', T_1, self.fluid_ref)
        cycle_states[1, 'P'], cycle_states[1, 'T'], cycle_states[1, 'H'], cycle_states[1, 'S'] = p_1, T_1, h_1, s_1

        # State 2: After superheater (superheated vapor, high pressure)
        T_1_sup = T_s_in - pp_1sup
        p_1_sup = p_1
        h_1_sup = PropsSI('H', 'P', p_1_sup, 'T', T_1_sup, self.fluid_ref)
        s_1_sup = PropsSI('S', 'P', p_1_sup, 'T', T_1_sup, self.fluid_ref)
        cycle_states[2, 'P'], cycle_states[2, 'T'], cycle_states[2, 'H'], cycle_states[2, 'S'] = p_1_sup, T_1_sup, h_1_sup, s_1_sup

        # State 3: After turbine (superheated vapor, low pressure)
        T_23 = T_cs_mid1 - pp_cs_mid1
        p_23 = PropsSI('P', 'Q', 1.0, 'T', T_23,self.fluid_ref)
        p_2 = p_23
        h_2_isentropic = PropsSI("H", "P", p_2, "S", s_1_sup, self.fluid_ref)
        h_2 = h_1_sup - self.eta_turbine * (h_1_sup - h_2_isentropic)
        T_2 = PropsSI('T', 'P', p_2, 'H', h_2, self.fluid_ref)
        s_2 = PropsSI('S', 'P', p_2, 'H', h_2, self.fluid_ref)
        cycle_states[3, 'P'], cycle_states[3, 'T'], cycle_states[3, 'H'], cycle_states[3, 'S'] = p_2, T_2, h_2, s_2
        # State 4: After desuperheater (saturated vapor, low pressure)
        h_23 = PropsSI('H', 'Q', 1.0, 'T', T_23, self.fluid_ref)
        s_23 = PropsSI('S', 'Q', 1.0, 'T', T_23, self.fluid_ref)
        cycle_states[4, 'P'], cycle_states[4, 'T'], cycle_states[4, 'H'], cycle_states[4, 'S'] = p_23, T_23, h_23, s_23
        # State 5: After condenser (saturated liquid, low pressure)
        T_3 = T_23
        p_3 = PropsSI('P', 'Q', 0.0, 'T', T_3,self.fluid_ref)
        h_3 = PropsSI('H', 'Q', 0.0, 'T', T_3, self.fluid_ref)
        s_3 = PropsSI('S', 'Q', 0.0, 'T', T_3, self.fluid_ref)
        cycle_states[5, 'P'], cycle_states[5, 'T'], cycle_states[5, 'H'], cycle_states[5, 'S'] = p_3, T_3, h_3, s_3
        # State 6: After pump (compressed liquid, high pressure)
        p_4 = p_0 
        h_4_isentropic = PropsSI("H", "P", p_4, "S", s_3, self.fluid_ref)
        h_4 = h_3 - self.eta_pump * (h_3 - h_4_isentropic)
        T_4 = PropsSI('T', 'P', p_4, 'H', h_4, self.fluid_ref)
        s_4 = PropsSI('S', 'P', p_4, 'H', h_4, self.fluid_ref)
        cycle_states[6, 'P'], cycle_states[6, 'T'], cycle_states[6, 'H'], cycle_states[6, 'S'] = p_4, T_4, h_4, s_4


        self.cycle_states = cycle_states
        self.fill_states()

    def extract_ts_data(self, states):
        """
        Extract T-s data from an external StateContainer.

        Parameters:
        -----------
        states : StateContainer
            The external state container to extract data from.

        Returns:
        --------
        np.ndarray
            A 2D array where each row contains [T, S].
        """
        # Extract T and S for all states
        temperatures = np.array([states[i, 'T'] for i in range(len(states))])  # Temperature in K
        entropies = np.array([states[i, 'S'] / 1000 for i in range(len(states))])  # Entropy in kJ/kg·K

        # Combine into a single array for easy processing
        ts_data = np.column_stack((temperatures, entropies))

        return ts_data

    def extract_ph_data(self, states):
        """
        Extract p-H data from an external StateContainer.

        Parameters:
        -----------
        states : StateContainer
            The external state container to extract data from.

        Returns:
        --------
        np.ndarray
            A 2D array where each row contains [P, H].
        """
        # Extract P and H for all states
        pressures = np.array([states[i, 'P'] / 1e6 for i in range(len(states))])  # Convert P to MPa
        enthalpies = np.array([states[i, 'H'] / 1000 for i in range(len(states))])  # Convert H to kJ/kg

        # Combine into a single array for easy processing
        ph_data = np.column_stack((pressures, enthalpies))

        return ph_data

    def eta_thermal(self):
        """
        Calculate the thermal efficiency of the cycle.

        Returns:
        --------
        float
            Thermal efficiency.
        """
        w_net = self.cycle_states[2].H - self.cycle_states[3].H - (self.cycle_states[6].H - self.cycle_states[5].H)
        q_boiler = self.cycle_states[2].H - self.cycle_states[0].H
        return w_net / q_boiler

    def exergy_efficiency(self, T_ref):
        """
        Calculate the exergy efficiency of the Rankine cycle.

        Parameters:
        -----------
        T_ref : float
            Reference temperature (K), usually the ambient temperature.

        Returns:
        --------
        float
            Exergy efficiency of the cycle.
        """
        # Exergy input (heat source)
        h_0 = self.cycle_states[0].H
        s_0 = self.cycle_states[0].S
        h_2 = self.cycle_states[2].H
        s_2 = self.cycle_states[2].S

        exergy_input = (h_2 - h_0) - T_ref * (s_2 - s_0)

        # Work output (net work)
        h_3 = self.cycle_states[3].H  # Turbine outlet
        h_6 = self.cycle_states[6].H  # Pump outlet
        h_5 = self.cycle_states[5].H  # Pump inlet

        turbine_work = h_2 - h_3
        pump_work = h_6 - h_5
        net_work = turbine_work - pump_work

        # Exergy efficiency
        exergy_eff = net_work / exergy_input
        return exergy_eff
    def hs_temp_out_net_power(self,p_hs_in, T_s_mid1, T_hs_in):
        """
        Calculate the temperature of the heat source outlet.

        Parameters:
        -----------
        m_dot_hs : float
            Mass flow rate of the heat source (kg/s).
        p_hs_in : float
            Pressure of the heat source inlet (Pa).

        Returns:
        --------
        float
            Temperature of the heat source outlet (K).
        """
        h_hs_mid = PropsSI('H', 'P', p_hs_in, 'T', T_s_mid1, self.fluid_ref)
        h_hs_in = PropsSI('H', 'P', p_hs_in, 'T', T_hs_in, self.fluid_ref)
        m_dot_wf = self.m_dot_hs * (h_hs_in - h_hs_mid) / (self.cycle_states[2].H - self.cycle_states[0].H)
        q_preheater_wf = m_dot_wf * (self.cycle_states[0].H - self.cycle_states[4].H)
        deltaH_preheater = q_preheater_wf / self.m_dot_hs
        h_hs_out = h_hs_mid - deltaH_preheater 
        T_hs_out = PropsSI('T', 'P', p_hs_in, 'H', h_hs_out, self.fluid_ref)
        

        w_net = m_dot_wf * (self.cycle_states[2].H - self.cycle_states[3].H - (self.cycle_states[6].H - self.cycle_states[5].H))
        return w_net, T_hs_out, m_dot_wf
    

        