from CoolProp.CoolProp import PropsSI
import CoolProp
from CoolProp.Plots.SimpleCycles import BaseCycle
import numpy as np
from StateContainer import StateContainer, StatePoint

class ORCSuperheat(BaseCycle):
    """
    Computes and stores thermodynamic states for a Rankine cycle with 
    superheating, a desuperheater, and a pump.
    """
    
    STATECOUNT = 7  # Number of states in the cycle
    STATECHANGE = [
        lambda inp: BaseCycle.state_change(inp, 'H', 'P', 0, ty1='lin', ty2='lin'),  # Isobaric heating (Heater)
        lambda inp: BaseCycle.state_change(inp, 'H', 'P', 1, ty1='lin', ty2='lin'),  # Isobaric heating (Superheater)
        lambda inp: BaseCycle.state_change(inp, 'S', 'P', 2, ty1='log', ty2='log'),  # Isentropic expansion (Turbine)
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
        m_dot_hs : float
            Mass flow rate of the heat source.
        graph_type : str, optional
            Type of thermodynamic graph to generate (default is "TS").
        """
        self.eta_turbine = eta_turbine
        self.eta_pump = eta_pump
        self.fluid_ref = fluid  
        self.m_dot_hs = m_dot_hs
        super().__init__(fluid, graph_type, **kwargs)

    def simple_solve(self, T_hs_in, T_pinch_hs_in, T_estimation_evap, 
                     T_pinch_hs_mid1, T_pinch_cs_mid1, T_cs_in, T_estimation_cndsr, SI=True):
        """
        Solve the thermodynamic states for the Rankine cycle with superheating, 
        desuperheating, and a pump.

        Parameters:
        -----------
        T_hs_in : float
            Heat source temperature (K).
        T_pinch_hs_in : float
            Superheater pinch point (K).
        T_estimation_evap : float
            Temperature difference for evaporator estimation (K). Must be > 0.
        T_pinch_hs_mid1 : float
            Middle source pinch point (K).
        T_pinch_cs_mid1 : float
            Middle cooling source pinch point (K).
        T_cs_in : float
            Inlet cooling source temperature (K).
        T_estimation_cndsr : float
            Temperature difference for condenser estimation (K). Must be > 0.
        SI : bool, optional
            Use SI units (default is True).

        Returns:
        --------
        StateContainer
            Container storing computed thermodynamic states.
        """
        cycle_states = StateContainer(unit_system=self._system)

        # Define temperature levels
        T_high = T_hs_in - T_pinch_hs_in
        T_hs_mid1 = T_hs_in - T_estimation_evap
        T_cs_mid1 = T_cs_in + T_estimation_cndsr
        T_low = T_cs_mid1 + T_pinch_cs_mid1

        # State 0: After preheater (saturated liquid, high pressure)
        T0, x0 = T_hs_mid1 - T_pinch_hs_mid1, 0
        p0 = PropsSI('P', 'Q', x0, 'T', T0, self.fluid_ref)
        h0 = PropsSI('H', 'Q', x0, 'T', T0, self.fluid_ref)
        s0 = PropsSI('S', 'Q', x0, 'T', T0, self.fluid_ref)
        cycle_states[0] = StatePoint(position="after preheater")
        cycle_states[0, 'P'], cycle_states[0, 'T'], cycle_states[0, 'H'], cycle_states[0, 'S'] = p0, T0, h0, s0
        

        # State 1: After heater (saturated vapor, high pressure)
        T1, x1 = T0, 1.0
        p1 = PropsSI('P', 'Q', x1, 'T', T1, self.fluid_ref)
        h1re = PropsSI('H', 'Q', x1, 'T', T1, self.fluid_ref)
        s1re = PropsSI('S', 'Q', x1, 'T', T1, self.fluid_ref)
        cycle_states[1] = StatePoint(position="after heater")
        cycle_states[1, 'P'], cycle_states[1, 'T'], cycle_states[1, 'H'], cycle_states[1, 'S'] = p1, T1, h1re, s1re

        # State 1sup: After superheater (superheated vapor, high pressure)
        T1sup, p1sup = T_high, p1
        h1sup = PropsSI('H', 'P', p1sup, 'T', T1sup, self.fluid_ref)
        s1sup = PropsSI('S', 'P', p1sup, 'T', T1sup, self.fluid_ref)
        cycle_states[2] = StatePoint(position="after superheater")
        cycle_states[2, 'P'], cycle_states[2, 'T'], cycle_states[2, 'H'], cycle_states[2, 'S'] = p1sup, T1sup, h1sup, s1sup

        # State 3: After condenser (saturated liquid, low pressure)
        T3, x3 = T_low, 0
        p3 = PropsSI('P', 'Q', x3, 'T', T3, self.fluid_ref)
        h3 = PropsSI('H', 'Q', x3, 'T', T3, self.fluid_ref)
        s3 = PropsSI('S', 'Q', x3, 'T', T3, self.fluid_ref)
        cycle_states[5] = StatePoint(position="after condenser")
        cycle_states[5, 'P'], cycle_states[5, 'T'], cycle_states[5, 'H'], cycle_states[5, 'S'] = p3, T3, h3, s3

        # State 2: After turbine (superheated vapor, low pressure)
        T23, x23 = T_low, 1.0
        p_23 = PropsSI('P', 'Q', x23, 'T', T23, self.fluid_ref)
        p2, s2is = p_23, s1sup
        h2is = PropsSI("H", "P", p2, "S", s2is, self.fluid_ref)
        h2re = h1sup - self.eta_turbine * (h1sup - h2is)
        T2 = PropsSI('T', 'P', p2, 'H', h2re, self.fluid_ref)
        s2re = PropsSI('S', 'P', p2, 'H', h2re, self.fluid_ref)
        cycle_states[3] = StatePoint(position="after turbine")
        cycle_states[3, 'P'], cycle_states[3, 'T'], cycle_states[3, 'H'], cycle_states[3, 'S'] = p2, T2, h2re, s2re

        # State 23: After desuperheater (saturated vapor, low pressure)
        h23re = PropsSI('H', 'Q', 1.0, 'T', T23, self.fluid_ref)
        s23re = PropsSI('S', 'Q', 1.0, 'T', T23, self.fluid_ref)
        cycle_states[4] = StatePoint(position="after desuperheater")
        cycle_states[4, 'P'], cycle_states[4, 'T'], cycle_states[4, 'H'], cycle_states[4, 'S'] = p_23, T23, h23re, s23re
        # State 6: After pump (compressed liquid, high pressure)
        p4 = p0 
        s4is = s3
        h4is = PropsSI("H", "P", p4, "S", s4is, self.fluid_ref)
        h4re = h3 - self.eta_pump * (h3 - h4is)
        T4 = PropsSI('T', 'P', p4, 'H', h4re, self.fluid_ref)
        s4re = PropsSI('S', 'P', p4, 'H', h4re, self.fluid_ref)
        cycle_states[6] = StatePoint(position="after pump")
        cycle_states[6, 'P'], cycle_states[6, 'T'], cycle_states[6, 'H'], cycle_states[6, 'S'] = p4, T4, h4re, s4re
        self.cycle_states = cycle_states
        self.fill_states()

        return cycle_states


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
        
        try:
            # **1. Compute Net Work Output (w_net)**
            w_net = self.cycle_states[2].H - self.cycle_states[3].H - (
                self.cycle_states[6].H - self.cycle_states[5].H
            )

            # **2. Compute Heat Input from Boiler (q_boiler)**
            q_boiler = max(self.cycle_states[2].H - self.cycle_states[0].H, 1e-6)  # Avoid division by zero

            # **3. Compute Thermal Efficiency**
            eta_th = w_net / q_boiler

        except KeyError as e:
            raise ValueError(f"Missing cycle state data: {e}")
        except Exception as e:
            raise ValueError(f"Error calculating thermal efficiency: {e}")

        return eta_th


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

        # **1. Validate Input Values**
        if T_ref <= 0:
            raise ValueError("Reference temperature (T_ref) must be greater than zero.")

        try:
            # **2. Compute Exergy Input (heat source)**
            h0re = self.cycle_states[0].H
            s0re = self.cycle_states[0].S
            h2re = self.cycle_states[2].H
            s2re = self.cycle_states[2].S

            exergy_input = max((h2re - h0re) - T_ref * (s2re - s0re), 1e-6)  # Avoid division by zero

            # **3. Compute Work Output (Net Work)**
            h3re = self.cycle_states[3].H  # Turbine outlet
            h_6 = self.cycle_states[6].H  # Pump outlet
            h_5 = self.cycle_states[5].H  # Pump inlet

            turbine_work = h2re - h3re
            pump_work = h_6 - h_5
            net_work = turbine_work - pump_work

            # **4. Compute Exergy Efficiency**
            exergy_eff = net_work / exergy_input

        except KeyError as e:
            raise ValueError(f"Missing cycle state data: {e}")
        except Exception as e:
            raise ValueError(f"Error calculating exergy efficiency: {e}")

        return exergy_eff

    def analyze_system(self, p_hs_in, T_estimation_evap, T_estimation_cndsr, T_hs_in, T_cs_in, p_cs_in):
        """
        Analyze the system by calculating:
        - Net work output.
        - Heat source outlet temperature.
        - Required mass flow rate of the working fluid.
        - Required mass flow rate of the cooling source.

        Parameters:
        -----------
        p_hs_in : float
            Pressure of the heat source inlet (Pa).
        T_estimation_evap : float
            Temperature difference for evaporator estimation (K). Must be > 0.
        T_estimation_cndsr : float
            Temperature difference for condenser estimation (K). Must be > 0.
        T_hs_in : float
            Temperature of the heat source inlet (K).
        T_cs_in : float
            Temperature of the cooling source inlet (K).
        p_cs_in : float
            Pressure of the cooling source inlet (Pa).

        Returns:
        --------
        tuple:
            - w_net (Net work output, W).
            - T_hs_out (Heat source outlet temperature, K).
            - m_dot_wf (Mass flow rate of working fluid, kg/s).
            - m_dot_cs (Mass flow rate of cooling source, kg/s).
        """

        # **1. Validate Input Values**
        if any(val <= 0 for val in [p_hs_in, T_hs_in, p_cs_in, T_cs_in]):
            raise ValueError("Pressure and temperature values must be greater than zero.")

        # **2. Compute Working Fluid Mass Flow Rate (m_dot_wf)**
        try:
            T_hs_mid1 = T_hs_in - T_estimation_evap
            h_hs_mid = PropsSI('H', 'P', p_hs_in, 'T', T_hs_mid1, self.fluid_ref)
            h_hs_in = PropsSI('H', 'P', p_hs_in, 'T', T_hs_in, self.fluid_ref)
            h_wf_in = self.cycle_states[0].H
            h_wf_out = self.cycle_states[2].H

            deltaH_wf = max(h_wf_out - h_wf_in, 1e-6)  # Avoid division by zero
            m_dot_wf = self.m_dot_hs * max(h_hs_in - h_hs_mid, 1e-6) / deltaH_wf

        except Exception as e:
            raise ValueError(f"Error calculating working fluid properties: {e}")

        # **3. Compute Heat Source Outlet Temperature (T_hs_out)**
        try:
            q_preheater_wf = m_dot_wf * (self.cycle_states[0].H - self.cycle_states[4].H)
            deltaH_preheater = q_preheater_wf / max(self.m_dot_hs, 1e-6)  # Avoid division by zero
            h_hs_out = h_hs_mid - deltaH_preheater
            T_hs_out = PropsSI('T', 'P', p_hs_in, 'H', h_hs_out, self.fluid_ref)

        except Exception as e:
            raise ValueError(f"Error calculating heat source outlet temperature: {e}")

        # **4. Compute Required Cooling Source Mass Flow Rate (m_dot_cs)**
        try:
            T_cs_mid1 = T_cs_in + T_estimation_cndsr
            h_cs_in = PropsSI('H', 'P', p_cs_in, 'T', T_cs_in, self.fluid_ref)
            h_cs_mid1 = PropsSI('H', 'P', p_cs_in, 'T', T_cs_mid1, self.fluid_ref)

            deltaH_23to3 = self.cycle_states[4].H - self.cycle_states[5].H
            deltaH_csintocsmid1 = max(h_cs_mid1 - h_cs_in, 1e-6)  # Avoid division by zero
            m_dot_cs = m_dot_wf * deltaH_23to3 / deltaH_csintocsmid1

        except Exception as e:
            raise ValueError(f"Error calculating cooling source mass flow rate: {e}")

        # **5. Compute Net Work Output (w_net)**
        try:
            w_net = m_dot_wf * (
                self.cycle_states[2].H - self.cycle_states[3].H - 
                (self.cycle_states[6].H - self.cycle_states[5].H)
            )
        except Exception as e:
            raise ValueError(f"Error calculating net work output: {e}")

        return w_net, T_hs_out, m_dot_wf, m_dot_cs
