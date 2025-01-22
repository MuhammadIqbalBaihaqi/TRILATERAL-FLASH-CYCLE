
import CoolProp
from CoolProp.CoolProp import PropsSI
from numpy import *
from matplotlib.pyplot import *
from TSdiagram import TSdiagram
from pHdiagram import pHdiagram
from CoolProp.Plots.SimpleCycles import StateContainer, BaseCycle
from ORCSuper import ORCSuperheat
from CoolProp.Plots import PropertyPlot
import warnings

T_s_in = 150 #dummies
pp_1sup = 5 #dummies
T_s_mid1 = 140  #dummies
pp_0 = 5 #dummies
T_cs_mid1 = 120 #dummies
T_cs_in = 273 #dummies
pp_cs_mid1 = 5 #dummies
eta_turbine = 0.5 #dummies
eta_pump = 0.5 #dummies
m_dot_hs = 1 #dummies
fluid = 'Methane'

rankine_cycle = ORCSuperheat(fluid, eta_turbine, eta_pump, m_dot_hs)
rankine_cycle.simple_solve(T_s_in, pp_1sup, T_s_mid1, pp_0, T_cs_mid1, pp_cs_mid1, T_cs_in)
warnings.filterwarnings('ignore')
rankine_cycle.steps = 20
halo = rankine_cycle.get_state_changes()

ts_data, ph_data = rankine_cycle.extract_ts_data(halo), rankine_cycle.extract_ph_data(halo)

plot_TSdiagram = TSdiagram(fluid, extracted_ts=ts_data)
plot_pHdiagram = pHdiagram(fluid, extracted_ph=ph_data)

show()
