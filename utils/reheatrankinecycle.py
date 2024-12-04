from CoolProp.CoolProp import PropsSI
from CoolProp.Plots.SimpleCycles import StateContainer, BaseCycle
""" Example usage:
from CoolProp.Plots import PropertyPlot
import warnings
warnings.filterwarnings('ignore')

# Define the system variables
fluid = 'Water'
p_high_in = 8e6          # High-pressure turbine inlet pressure (Pa)
T_high_in = 480 + 273.15 # High-pressure turbine inlet temperature (K)
p_high_out = 7e5         # High-pressure turbine outlet pressure (Pa)
T_reheat_in = 440 + 273.15 # Reheat temperature before low-pressure turbine (K)
eta_turbine_high = 1.0   # High-pressure turbine efficiency
eta_turbine_low = 1.0    # Low-pressure turbine efficiency
p_cond = 0.008e6         # Condenser pressure (Pa)
q_cond = 0.0             # Quality at condenser outlet (saturated liquid)
eta_pump = 1.0           # Pump efficiency

# Initialize the Rankine cycle with reheating
rankine_cycle = ReheatRankineCycle(fluid, eta_turbine_high, eta_turbine_low, eta_pump)

# Compute the thermodynamic cycle
container = rankine_cycle.simple_solve(
    p_high_in=p_high_in,
    T_high_in=T_high_in,
    p_high_out=p_high_out,
    T_reheat_in=T_reheat_in,
    p_cond=p_cond,
    q_cond=q_cond
)

# Retrieve and display the states
rankine_cycle.steps = 200
halo = rankine_cycle.get_state_changes()

plot = PropertyPlot('Water', 'TS')
plot.calc_isolines(CoolProp.iP, num=15)
plot.calc_isolines(CoolProp.iQ, num=15)
plot.draw_process(halo)
plt.close(rankine_cycle.figure)
plt.ioff()

# Show the plot
plot.show()
"""
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

    def __init__(self, fluid, eta_turbine_high, eta_turbine_low, eta_pump, graph_type="TS", **kwargs):
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
        w_net = self.cycle_states[2].H - self.cycle_states[3].H + self.cycle_states[4].H - self.cycle_states[5].H- (self.cycle_states[1].H - self.cycle_states[0].H)
        q_boiler = self.cycle_states[2].H - self.cycle_states[1].H + self.cycle_states[4].H - self.cycle_states[3].H
        return w_net / q_boiler
