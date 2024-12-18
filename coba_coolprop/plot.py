import CoolProp
from CoolProp.CoolProp import PropsSI
from CoolProp.Plots import PropertyPlot
from CoolProp.Plots import SimpleRankineCycle
import matplotlib.pyplot as plt
import warnings


class RankineCyclePlotter:
    def __init__(self, fluid, diagram_type='TS', unit_system='SI'):
        """
        Initialize the RankineCyclePlotter class.

        Args:
            fluid (str): The fluid name in CoolProp format.
            diagram_type (str): The type of plot (e.g., 'TS', 'PH').
            unit_system (str): The unit system (default is 'SI').
        """
        self.fluid = fluid
        self.diagram_type = diagram_type
        self.unit_system = unit_system
        self.pp = PropertyPlot(f'HEOS::{fluid}', diagram_type, unit_system=unit_system)

    def add_isolines(self, variable, num=10):
        """
        Add isolines to the property plot.

        Args:
            variable (str): The property for isolines (e.g., CoolProp.iQ, CoolProp.iP).
            num (int): Number of isolines to calculate.
        """
        self.pp.calc_isolines(variable, num=num)

    def setup_rankine_cycle(self, T0, T2, p0, p2, eff_pump, eff_turbine):
        """
        Set up and solve the Simple Rankine Cycle.

        Args:
            T0 (float): Condenser temperature in Kelvin.
            T2 (float): Boiler temperature in Kelvin.
            eff_pump (float): Isentropic efficiency of the pump.
            eff_turbine (float): Isentropic efficiency of the turbine.

        Returns:
            tuple: A tuple containing the Rankine cycle object and state changes.
        """
        cycle = SimpleRankineCycle(f'HEOS::{self.fluid}', self.diagram_type, unit_system=self.unit_system)

        # Calculate pressures at condenser and boiler states
        state = cycle.state
        state.update(CoolProp.QT_INPUTS, 0.0, T0 - 10)
        state.update(CoolProp.QT_INPUTS, 1.0, T2 + 15)

        # Solve the cycle
        cycle.simple_solve(T0, p0, T2, p2, eff_pump, eff_turbine, SI=True)
        cycle.steps = 50
        sc = cycle.get_state_changes()
    
        return cycle, sc
    

    def calculate_efficiency(self, p0, p2, T0, T2, eff_pump, eff_turbine):
        """
        Calculate the thermal efficiency of the Rankine cycle using PropsSI.

        Args:
            p0 (float): Condenser pressure in Pa.
            p2 (float): Boiler pressure in Pa.
            T0 (float): Condenser temperature in Kelvin.
            T2 (float): Boiler temperature in Kelvin.
            eff_pump (float): Isentropic efficiency of the pump.
            eff_turbine (float): Isentropic efficiency of the turbine.

        Returns:
            float: The thermal efficiency of the cycle.
        """
        # State 1: Saturated liquid at condenser pressure
        h1 = PropsSI('H', 'P', p0, 'Q', 0, self.fluid)

        # State 2: Compressed liquid at boiler pressure (isentropic pump)
        s1 = PropsSI('S', 'P', p0, 'Q', 0, self.fluid)
        h2s = PropsSI('H', 'P', p2, 'S', s1, self.fluid)
        h2 = h1 + (h2s - h1) / eff_pump

        # State 3: Superheated vapor at boiler pressure and temperature
        h3 = PropsSI('H', 'P', p2, 'T', T2, self.fluid)

        # State 4: Expanded vapor at condenser pressure (isentropic turbine)
        s3 = PropsSI('S', 'P', p2, 'T', T2, self.fluid)
        h4s = PropsSI('H', 'P', p0, 'S', s3, self.fluid)
        h4 = h3 - eff_turbine * (h3 - h4s)

        # Work and heat calculations
        w_turbine = h3 - h4  # Turbine work
        w_pump = h2 - h1  # Pump work
        q_in = h3 - h2  # Heat input

        # Efficiency calculation
        efficiency = (w_turbine - w_pump) / q_in * 100  # Percentage

        return efficiency
    
    def plot_rankine_cycle(self, T0, T2, p0,p2, eff_pump, eff_turbine):
        """
        Plot the Rankine Cycle on a T-S diagram.

        Args:
            T0 (float): Condenser temperature in Kelvin.
            T2 (float): Boiler temperature in Kelvin.
            eff_pump (float): Isentropic efficiency of the pump.
            eff_turbine (float): Isentropic efficiency of the turbine.
        """
        # Add quality and pressure isolines
        self.add_isolines(CoolProp.iQ, num=11)
        self.add_isolines(CoolProp.iP, num=5)

        # Set up and solve the Rankine Cycle
        cycle, sc = self.setup_rankine_cycle(T0, T2, p0, p2, eff_pump, eff_turbine)

        # Draw the cycle on the T-S diagram
        self.pp.draw_process(sc)

        plt.close(cycle.figure)
        # Finalize the plot
        self.pp.show()
        
        # Calculate and print efficiency
        efficiency = self.calculate_efficiency(p0, p2, T0, T2, eff_pump, eff_turbine)
        print(f"Thermal Efficiency: {efficiency:.2f}%")

        warnings.filterwarnings("ignore", category=UserWarning)  # Replace UserWarning with the relevant warning type


# Example Usage
if __name__ == "__main__":
    fluid = 'R134a'
    p0 = 10000
    p2 = 1000000
    T2 = 350  # Boiler temperature in Kelvin
    T0 = PropsSI('T', 'Q', 1, 'P', p0, fluid) - .1
    eff_pump = 0.8  # Pump efficiency
    eff_turbine = 0.7  # Turbine efficiency

    plotter = RankineCyclePlotter(fluid)
    plotter.plot_rankine_cycle(T0, T2, p0, p2, eff_pump, eff_turbine)
