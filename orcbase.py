import CoolProp.CoolProp as CP
import warnings

warnings.filterwarnings("ignore")

class ORCBase():
    """Initialize the ORC system with the given parameters."""
    def __init__(self, eta_exp, eta_pmp, T_hs_in, T_cs_in, fluid_ref):
        self.fluid      = fluid_ref
        self.eta_exp    = eta_exp
        self.eta_pmp    = eta_pmp
        self.T_cs_in    = T_cs_in
        self.T_hs_in    = T_hs_in


    def simple_solve(self, T_estimation_cs, T_estimation_hs, T_pinch_hs, T_pinch_cs):
        """Solve the ORC cycle with ideal and real conditions.  """
        # Define cooling source temperature
        self.T_cs_mid1 = self.T_cs_in + T_estimation_cs
        self.T_hs_mid1 = self.T_hs_in - T_estimation_hs

        # Define high and low temperature
        T_high  = self.T_hs_mid1 - T_pinch_hs
        T_lower = self.T_cs_mid1 + T_pinch_cs

        ##########################################################
        #################### IDEAL CONDITIONS ####################
        ##########################################################
        # state 0 (Input: T0 dan X0)
        T0 = T_high
        X0 = 0
        P0 = CP.PropsSI('P', 'T', T0, 'Q', X0, self.fluid)
        H0 = CP.PropsSI('H', 'T', T0, 'Q', X0, self.fluid)
        S0 = CP.PropsSI('S', 'T', T0, 'Q', X0, self.fluid)

        # state 1 (Input: T1 dan X1)
        T1 = T0
        X1 = 1
        P1 = CP.PropsSI('P', 'T', T1, 'Q', X1, self.fluid)
        H1 = CP.PropsSI('H', 'T', T1, 'Q', X1, self.fluid)
        S1 = CP.PropsSI('S', 'T', T1, 'Q', X1, self.fluid)

        # state 3 (Input: T3 dan X3)
        T3 = T_lower
        X3 = 0
        P3 = CP.PropsSI('P', 'T', T3, 'Q', X3, self.fluid)
        H3 = CP.PropsSI('H', 'T', T3, 'Q', X3, self.fluid)
        S3 = CP.PropsSI('S', 'T', T3, 'Q', X3, self.fluid)

        # state 2 (Input: p2 dan s2)
        S2 = S1
        P2 = P3
        H2 = CP.PropsSI('H', 'S', S2, 'P', P2, self.fluid)
        T2 = CP.PropsSI('T', 'S', S2, 'P', P2, self.fluid)

        # state 23 (Input: X23 dan T2)
        P23 = P3
        X23 = 1
        T23 = T_lower
        H23 = CP.PropsSI('H', 'T', T23, 'Q', X23, self.fluid)
        S23 = CP.PropsSI('S', 'T', T23, 'Q', X23, self.fluid)

        # state 4 (Input: s4 dan p4)
        S4 = S3
        P4 = P0
        H4 = CP.PropsSI('H', 'S', S4, 'P', P4, self.fluid)
        T4 = CP.PropsSI('T', 'S', S4, 'P', P4, self.fluid)

        ##########################################################
        #################### REAL CONDITIONS #####################
        ##########################################################
        # state 0 (Input: T0 dan X0) before evaporator
        T0re = T0-1
        X0re = X0
        P0re = CP.PropsSI('P', 'T', T0re, 'Q', X0re, self.fluid)
        H0re = CP.PropsSI('H', 'T', T0re, 'Q', X0re, self.fluid)
        S0re = CP.PropsSI('S', 'T', T0re, 'Q', X0re, self.fluid)

        # state 1 (Input: T1 dan X1) before expander
        T1re = T0re-1
        X1re = X1
        P1re = CP.PropsSI('P', 'T', T1re, 'Q', X1re, self.fluid)
        H1re = CP.PropsSI('H', 'T', T1re, 'Q', X1re, self.fluid)
        S1re = CP.PropsSI('S', 'T', T1re, 'Q', X1re, self.fluid)

        # state 3 (Input: T3 dan X3) after condenser
        T3re = T3-2
        X3re = 0
        P3re = CP.PropsSI('P', 'T', T3re, 'Q', X3re, self.fluid)
        H3re = CP.PropsSI('H', 'T', T3re, 'Q', X3re, self.fluid)
        S3re = CP.PropsSI('S', 'T', T3re, 'Q', X3re, self.fluid)

        # state 2 (Input: p2 dan s2) after expander
        P2re    = P3 #perlu konfirmasi apakah referensi p2 real menggunakan p3 ideal atau p3 real. saya memakai 
                    # P2re = P3 karena P3 real sudah mengalami temperature drop 
                    # yang berarti juga adanya pressure drop. jadi, menurut saya p2 tdk sama dgn p3
                    # atau justru p2re = p23re 
        S2re_is = S1re
        H2re_is = CP.PropsSI('H', 'S', S2re_is, 'P', P2re, self.fluid)
        H2re    = H1re - self.eta_exp * (H1re - H2re_is)
        S2re    = CP.PropsSI('S', 'H', H2re, 'P', P2re, self.fluid)
        T2re    = CP.PropsSI('T', 'H', H2re, 'P', P2re, self.fluid)

        # state 23 (Input: X23 dan T2) after desuperheater
        X23re = 1
        T23re = T23 - 1
        P23re = CP.PropsSI('P', 'T', T23re, 'Q', X23re, self.fluid)
        H23re = CP.PropsSI('H', 'T', T23re, 'Q', X23re, self.fluid)
        S23re = CP.PropsSI('S', 'T', T23re, 'Q', X23re, self.fluid)

        # state 4 (Input: s4 dan p4) after pump
        S4re_is = S3re
        P4re    = P0re
        H4re_is = CP.PropsSI('H', 'S', S4re_is, 'P', P4re, self.fluid)
        H4re    = H3re + (H4re_is - H3re)/self.eta_pmp
        S4re    = CP.PropsSI('S', 'H', H4re, 'P', P4re, self.fluid)
        T4re    = CP.PropsSI('T', 'H', H4re, 'P', P4re, self.fluid)

        self.H_is_list = [H0, H1, H2, H23, H3, H4]
        self.s_is_list = [S0, S1, S2, S23, S3, S4]
        self.p_is_list = [P0, P1, P2, P23, P3, P4]
        self.T_is_list = [T0, T1, T2, T23, T3, T4]

        self.H_re_list = [H0re, H1re, H2re, H23re, H3re, H4re]
        self.s_re_list = [S0re, S1re, S2re, S23re, S3re, S4re]
        self.p_re_list = [P0re, P1re, P2re, P23re, P3re, P4re]
        self.T_re_list = [T0re, T1re, T2re, T23re, T3re, T4re]

    def calc_thermal_efficiency(self, real=True):
        """
        Calculate the thermal efficiency of the ORC cycle.
        If real=True, use real condition enthalpies, else use ideal.
        """
        if real:
            H1, H2, H4, H3 = self.H_re_list[1], self.H_re_list[2], self.H_re_list[5], self.H_re_list[4]
        else:
            H1, H2, H4, H3 = self.H_is_list[1], self.H_is_list[2], self.H_is_list[5], self.H_is_list[4]

        W_turbine = H1 - H2
        W_pump = H4 - H3
        Q_in = H1 - H4

        if Q_in == 0:
            return 0.0

        eta_thermal = (W_turbine - W_pump) / Q_in
        return eta_thermal

    # def calc_wf_mass_flow_from_heat_source(
    #     self, hs_mass_flow, p_hs, p_cs,
    #     hs_fluid_ref='Water', cs_fluid_ref='Water', real=True):
    #     """
    #     Calculate the ORC working fluid mass flow rate given heat source mass flow rate.
    #     Args:
    #         hs_mass_flow (float): Mass flow rate of heat source [kg/s]
    #         p_hs (float): Pressure of heat source [Pa]
    #         p_cs (float): Pressure of cooling source [Pa]
    #         hs_fluid_ref (str): Reference fluid for heat source (default: 'Water')
    #         cs_fluid_ref (str): Reference fluid for cooling source (default: 'Water')
    #         real (bool): Use real or ideal cycle properties
    #     Returns:
    #         tuple: (
    #             working fluid mass flow rate [kg/s],
    #             heat source outlet enthalpy [J/kg],
    #             heat source outlet temperature [K],
    #             heat source mid enthalpy [J/kg],
    #             heat source outlet entropy [J/kg·K],
    #             heat source mid entropy [J/kg·K],
    #             cooling source outlet enthalpy [J/kg],
    #             cooling source outlet temperature [K],
    #             cooling source outlet entropy [J/kg·K]
    #         )
    #     """
    #     # Heat source properties
    #     h_hs_in = CP.PropsSI('H', 'T', self.T_hs_in, 'P', p_hs, hs_fluid_ref)
    #     h_hs_mid1 = CP.PropsSI('H', 'T', self.T_hs_mid1, 'P', p_hs, hs_fluid_ref)
    #     s_hs_mid1 = CP.PropsSI('S', 'H', h_hs_mid1, 'P', p_hs, hs_fluid_ref)

    #     # Working fluid enthalpy (ideal or real)
    #     H = self.H_re_list if real else self.H_is_list

    #     # Heat transferred from heat source
    #     Q_give = hs_mass_flow * (h_hs_in - h_hs_mid1)
    #     Q_accept_spec = H[1] - H[0]
    #     wf_mass_flow = Q_give / Q_accept_spec if Q_accept_spec != 0 else 0.0

    #     # Heat source outlet enthalpy and properties
    #     h_hs_out = h_hs_mid1 - (wf_mass_flow * (H[0] - H[5]) / hs_mass_flow) if hs_mass_flow != 0 else h_hs_mid1
    #     T_hs_out = CP.PropsSI('T', 'H', h_hs_out, 'P', p_hs, hs_fluid_ref)
    #     s_hs_out = CP.PropsSI('S', 'H', h_hs_out, 'P', p_hs, hs_fluid_ref)

    #     # Cooling source outlet properties
    #     h_cs_mid1 = CP.PropsSI('H', 'T', self.T_cs_mid1, 'P', p_cs, cs_fluid_ref)
    #     h_cs_out = h_cs_mid1 + (wf_mass_flow * (H[2] - H[3]) / hs_mass_flow) if hs_mass_flow != 0 else h_cs_mid1
    #     T_cs_out = CP.PropsSI('T', 'H', h_cs_out, 'P', p_cs, cs_fluid_ref)
    #     s_cs_out = CP.PropsSI('S', 'H', h_cs_out, 'P', p_cs, cs_fluid_ref)

    #     return (
    #         wf_mass_flow, h_hs_out, T_hs_out, h_hs_mid1,
    #         s_hs_out, s_hs_mid1, h_cs_out, T_cs_out, s_cs_out
    #     )
    
    def calc_wf_mass_flow_from_heat_source(
        self, hs_mass_flow, p_hs, p_cs,
        hs_fluid_ref='Water', cs_fluid_ref='Water', real=True):
        """
        Calculate the ORC working fluid mass flow rate given heat source mass flow rate.

        Args:
            hs_mass_flow (float): Mass flow rate of heat source [kg/s]
            p_hs (float): Pressure of heat source [Pa]
            p_cs (float): Pressure of cooling source [Pa]
            hs_fluid_ref (str): Reference fluid for heat source (default: 'Water')
            cs_fluid_ref (str): Reference fluid for cooling source (default: 'Water')
            real (bool): Use real or ideal cycle properties

        Returns:
            dict: {
                "wf": {"mass_flow": ..., },
                "hs": {"h_mid1": ..., "h_out": ..., "T_out": ..., "s_mid1": ..., "s_out": ...},
                "cs": {"h_out": ..., "T_out": ..., "s_out": ...}
            }
        """
        # Heat source properties
        h_hs_in = CP.PropsSI('H', 'T', self.T_hs_in, 'P', p_hs, hs_fluid_ref)
        h_hs_mid1 = CP.PropsSI('H', 'T', self.T_hs_mid1, 'P', p_hs, hs_fluid_ref)
        s_hs_mid1 = CP.PropsSI('S', 'H', h_hs_mid1, 'P', p_hs, hs_fluid_ref)

        # Working fluid enthalpy (ideal or real)
        H = self.H_re_list if real else self.H_is_list

        # Heat transferred from heat source
        Q_give = hs_mass_flow * (h_hs_in - h_hs_mid1)
        Q_accept_spec = H[1] - H[0]
        wf_mass_flow = Q_give / Q_accept_spec #if Q_accept_spec != 0 else 0.0



        # Heat source outlet enthalpy and properties
        h_hs_out = h_hs_mid1 - (wf_mass_flow * (H[0] - H[5]) / hs_mass_flow) #if hs_mass_flow != 0 else h_hs_mid1
        T_hs_out = CP.PropsSI('T', 'H', h_hs_out, 'P', p_hs, hs_fluid_ref)
        s_hs_out = CP.PropsSI('S', 'H', h_hs_out, 'P', p_hs, hs_fluid_ref)

        # Cooling source outlet properties
        h_cs_in = CP.PropsSI('H', 'T', self.T_cs_in, 'P', p_cs, cs_fluid_ref)
        h_cs_mid1 = CP.PropsSI('H', 'T', self.T_cs_mid1, 'P', p_cs, cs_fluid_ref)
        # Calculate cooling source mass flow rate
        Q_give_cs = (H[3] - H[4]) * wf_mass_flow
        Q_accept_cs = h_cs_mid1 - h_cs_in
        cs_mass_flow = Q_give_cs / Q_accept_cs
        h_cs_out = h_cs_mid1 + (wf_mass_flow * (H[2] - H[3]) / cs_mass_flow) #if hs_mass_flow != 0 else h_cs_mid1
        T_cs_out = CP.PropsSI('T', 'H', h_cs_out, 'P', p_cs, cs_fluid_ref)
        s_cs_out = CP.PropsSI('S', 'H', h_cs_out, 'P', p_cs, cs_fluid_ref)



        return {
            "wf": {
                "mass_flow": wf_mass_flow
            },
            "hs": {
                "h_in": h_hs_in,
                "h_mid1": h_hs_mid1,
                "h_out": h_hs_out,
                "T_out": T_hs_out,
                "T_mid1": self.T_hs_mid1,
                "s_mid1": s_hs_mid1,
                "s_out": s_hs_out
            },
            "cs": {
                "mass_flow": cs_mass_flow,
                "h_in": h_cs_in,
                "h_mid1": h_cs_mid1,
                "h_out": h_cs_out,
                "T_mid1": self.T_cs_mid1,
                "T_out": T_cs_out,
                "s_out": s_cs_out
            }
        }

    # def calc_heat_source_properties(self, hs_mass_flow, p_hs, hs_fluid_ref='Water'):
    #     """
    #     Calculate heat source output properties after heat exchange.
    #     Args:
    #         hs_mass_flow (float): Mass flow rate of heat source [kg/s]
    #         p_hs (float): Pressure of heat source [Pa]
    #         hs_fluid_ref (str): Reference fluid for heat source (default: 'Water')
    #     Returns:
    #         dict: {'T_out': temperature output [K], 'h_out': enthalpy output [J/kg]}
    #     """
    #     T_out = self.T_hs_mid1
    #     h_out = CP.PropsSI('H', 'T', T_out, 'P', p_hs, hs_fluid_ref)
    #     return {'T_out': T_out, 'h_out': h_out}

    # def calc_specific_exergy_destruction(self, T0, wf_mass_flow, hs_mass_flow, cs_mass_flow, h_cs_in, h_cs_mid1, h_cs_out, h_hs_mid1, h_hs_out, s_hs_mid1, s_hs_out, h_hs_in, s_hs_in, s_cs_mid1, s_cs_out, s_cs_in, real=True):
    #     """
    #     Calculate specific exergy destruction for each component in the ORC cycle.
    #     Args:
    #         T0 (float): Ambient/reference temperature [K]
    #         wf_mass_flow (float): Working fluid mass flow rate [kg/s]
    #         hs_mass_flow (float): Heat source mass flow rate [kg/s]
    #         cs_mass_flow (float): Cooling source mass flow rate [kg/s]
    #         h_cs_in (float): Enthalpy of cooling source at inlet [J/kg]
    #         h_cs_mid1 (float): Enthalpy of cooling source at mid-point [J/kg]
    #         h_cs_out (float): Enthalpy of cooling source at outlet [J/kg]
    #         s_cs_mid1 (float): Entropy of cooling source at mid-point [J/kg·K]
    #         s_cs_out (float): Entropy of cooling source at outlet [J/kg·K]
    #         h_hs_mid1 (float): Enthalpy of heat source at mid-point [J/kg]
    #         h_hs_out (float): Enthalpy of heat source at outlet [J/kg]
    #         s_hs_mid1 (float): Entropy of heat source at mid-point [J/kg·K]
    #         s_hs_out (float): Entropy of heat source at outlet [J/kg·K]
    #         h_hs_in (float): Enthalpy of heat source at inlet [J/kg]
    #         s_hs_in (float): Entropy of heat source at inlet [J/kg·K]
    #         s_cs_in (float): Entropy of cooling source at inlet [J/kg·K]
    #         real (bool): Use real or ideal cycle properties
    #     Returns:
    #         dict: Specific exergy destruction for each component [J/kg]
    #     """
    #     if real:
    #         H = self.H_re_list
    #         S = self.s_re_list
    #     else:
    #         H = self.H_is_list
    #         S = self.s_is_list

    #     # Liquid Heater
    #     flow_ex_delta_hs_liq = hs_mass_flow * ((h_hs_mid1 - h_hs_out) - T0 * (s_hs_mid1 - s_hs_out))
    #     flow_ex_delta_wf_liq = wf_mass_flow * (H[5] - H[0] - T0 * (S[5] - S[0]))
    #     ex_dest_heater = flow_ex_delta_hs_liq + flow_ex_delta_wf_liq
    #     # Evaporator
    #     flow_ex_delta_hs_evap = hs_mass_flow * ((h_hs_in - h_hs_mid1) - T0 * (s_hs_in - s_hs_mid1))
    #     flow_ex_delta_wf_evap = wf_mass_flow * (H[0] - H[1] - T0 * (S[0] - S[1]))
    #     ex_dest_evaporator = flow_ex_delta_hs_evap + flow_ex_delta_wf_evap
    #     # Turbine
    #     ex_dest_turbine = (- T0 * (S[1] - S[2])) * wf_mass_flow
    #     # Desuperheater
    #     flow_ex_delta_cs_desup = cs_mass_flow * (h_cs_mid1 - h_cs_out - T0 * (s_cs_mid1 - s_cs_out))
    #     flow_ex_delta_wf_desup = wf_mass_flow * (H[2] - H[3] - T0 * (S[2] - S[3]))
    #     ex_dest_desuperheater = flow_ex_delta_cs_desup + flow_ex_delta_wf_desup
    #     # Condenser
    #     flow_ex_delta_cs_cond = cs_mass_flow * ((h_cs_in - h_cs_mid1) - T0 * (s_cs_in - s_cs_mid1))
    #     flow_ex_delta_wf_cond = wf_mass_flow * (H[3] - H[4] - T0 * (S[3] - S[4]))
    #     ex_dest_condenser = flow_ex_delta_cs_cond + flow_ex_delta_wf_cond
    #     # Pump
    #     ex_dest_pump = (- T0 * (S[4] - S[5])) * wf_mass_flow

    #     return {
    #         'liquid heater': ex_dest_heater,
    #         'evaporator': ex_dest_evaporator,
    #         'turbine': ex_dest_turbine,
    #         'desuperheater': ex_dest_desuperheater,
    #         'condenser': ex_dest_condenser,
    #         'pump': ex_dest_pump
    #     }

    def calc_specific_exergy_destruction(self, T0, wf_mass_flow, hs_mass_flow, cs_mass_flow, hs, cs, real=True):
        """
        Args:
            T0 (float): Dead state temperature [K]
            wf_mass_flow (float): Working fluid mass flow rate [kg/s]
            hs_mass_flow (float): Heat source mass flow rate [kg/s]
            cs_mass_flow (float): Cooling source mass flow rate [kg/s]
            hs (dict): Heat source enthalpy & entropy {"h": [h_in, h_mid1, h_out], "s": [s_in, s_mid1, s_out]}
            cs (dict): Cooling source enthalpy & entropy {"h": [h_in, h_mid1, h_out], "s": [s_in, s_mid1, s_out]}
            real (bool): Use real or ideal cycle properties
        """

        if real:
            H = self.H_re_list
            S = self.s_re_list
        else:
            H = self.H_is_list
            S = self.s_is_list

        # Unpack for readability
        h_hs_in, h_hs_mid1, h_hs_out = hs["h"]
        s_hs_in, s_hs_mid1, s_hs_out = hs["s"]

        h_cs_in, h_cs_mid1, h_cs_out = cs["h"]
        s_cs_in, s_cs_mid1, s_cs_out = cs["s"]

        # Liquid Heater
        flow_ex_delta_hs_liq = hs_mass_flow * ((h_hs_mid1 - h_hs_out) - T0 * (s_hs_mid1 - s_hs_out))
        flow_ex_delta_wf_liq = wf_mass_flow * (H[5] - H[0] - T0 * (S[5] - S[0]))
        ex_dest_heater = flow_ex_delta_hs_liq + flow_ex_delta_wf_liq

        # Evaporator
        flow_ex_delta_hs_evap = hs_mass_flow * ((h_hs_in - h_hs_mid1) - T0 * (s_hs_in - s_hs_mid1))
        flow_ex_delta_wf_evap = wf_mass_flow * (H[0] - H[1] - T0 * (S[0] - S[1]))
        ex_dest_evaporator = flow_ex_delta_hs_evap + flow_ex_delta_wf_evap

        # Turbine
        ex_dest_turbine = (- T0 * (S[1] - S[2])) * wf_mass_flow

        # Desuperheater
        flow_ex_delta_cs_desup = cs_mass_flow * (h_cs_mid1 - h_cs_out - T0 * (s_cs_mid1 - s_cs_out))
        flow_ex_delta_wf_desup = wf_mass_flow * (H[2] - H[3] - T0 * (S[2] - S[3]))
        ex_dest_desuperheater = flow_ex_delta_cs_desup + flow_ex_delta_wf_desup

        # Condenser
        flow_ex_delta_cs_cond = cs_mass_flow * ((h_cs_in - h_cs_mid1) - T0 * (s_cs_in - s_cs_mid1))
        flow_ex_delta_wf_cond = wf_mass_flow * (H[3] - H[4] - T0 * (S[3] - S[4]))
        ex_dest_condenser = flow_ex_delta_cs_cond + flow_ex_delta_wf_cond

        # Pump
        ex_dest_pump = (- T0 * (S[4] - S[5])) * wf_mass_flow

        return {
            'liquid heater': ex_dest_heater,
            'evaporator': ex_dest_evaporator,
            'turbine': ex_dest_turbine,
            'desuperheater': ex_dest_desuperheater,
            'condenser': ex_dest_condenser,
            'pump': ex_dest_pump
        }
