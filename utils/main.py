import CoolProp
from CoolProp.CoolProp import PropsSI
from numpy import *
from matplotlib.pyplot import *
from TSdiagram import TSdiagram
from pHdiagram import pHdiagram
from CoolProp.Plots.SimpleCycles import BaseCycle
from ORCsuper_rev import ORCSuperheat
from CoolProp.Plots import PropertyPlot
import warnings


T_hs_in = 100 + 273.15 #dummies
T_estimation_evap = 20
T_pinch_hs_mid1, T_pinch_hs_in = 5,5
T_estimation_cndsr = 30 #dummies
T_cs_in = -15+273.15 #dummies
T_pinch_cs_mid1 = 5 #dummies
eta_turbine = 0.5 #dummies
eta_pump = 0.5 #dummies
m_dot_hs = 1 #dummies
fluid = 'DME' #T_hs_in, T_pinch_hs_in, T_estimation_evap, T_pinch_hs_mid1, T_pinch_cs_mid1, T_cs_in, T_estimation_cndsr

rankine_cycle = ORCSuperheat(fluid, eta_turbine, eta_pump, m_dot_hs)
close()
states = rankine_cycle.simple_solve(T_hs_in, T_pinch_hs_in, T_estimation_evap, T_pinch_hs_mid1, T_pinch_cs_mid1, T_cs_in, T_estimation_cndsr)
warnings.filterwarnings('ignore')
# rankine_cycle.steps = 1
# halo = rankine_cycle.get_state_changes()


# ts_data, ph_data = rankine_cycle.extract_ts_data(halo), rankine_cycle.extract_ph_data(halo)

# plot_TSdiagram = TSdiagram(fluid, extracted_ts=ts_data)
# plot_pHdiagram = pHdiagram(fluid, extracted_ph=ph_data)
print(states)
# show()
