from orcbase import ORCBase
from TSdiagram import TSdiagram
import matplotlib.pyplot as plt
import numpy as np
import exergy_method as exer

# def main():
#     # bikin object ORC
#     orc = ORCBase(
#         eta_exp = 0.90,             # efisiensi expander
#         eta_pmp = 0.90,              # efisiensi pump
#         T_hs_in = 100 + 273.15,     # hot source inlet, misal 450°C (K)
#         T_cs_in = 10 + 273.15,      # cold source inlet, misal 25°C (K)
#         fluid_ref = "Water"        # working fluid
#     )

#     # solve ORC
#     orc.simple_solve(
#         T_estimation_cs = 5,
#         T_estimation_hs = 5,
#         T_pinch_hs = 5,
#         T_pinch_cs = 5
#     )

#     # print hasil
#     print("=== IDEAL CYCLE RESULTS ===")
#     print("Temperatures (K):", orc.T_is_list)
#     print("Pressures (Pa):", orc.p_is_list)
#     print("Entropies (J/kg-K):", orc.s_is_list)
#     print("Enthalpies (J/kg):", orc.H_is_list)

#     print("\n=== REAL CYCLE RESULTS ===")
#     print("Temperatures (K):", orc.T_re_list)
#     print("Pressures (Pa):", orc.p_re_list)
#     print("Entropies (J/kg-K):", orc.s_re_list)
#     print("Enthalpies (J/kg):", orc.H_re_list)

#     # hitung efisiensi termal ideal dan nyata
#     eta_thermal_ideal = orc.calc_thermal_efficiency(real=False)
#     eta_thermal_real = orc.calc_thermal_efficiency(real=True)

#     print("\n=== THERMAL EFFICIENCY ===")
#     print(f"Ideal cycle efficiency: {eta_thermal_ideal:.4f}")
#     print(f"Real cycle efficiency: {eta_thermal_real:.4f}")
# Dummy data for demonstration
orc = ORCBase(
    eta_exp=0.6,
    eta_pmp=0.7,
    T_hs_in=122 + 273.15,    # K
    T_cs_in=10 + 273.15,    # K
    fluid_ref='R245fa'
)

# Solve the cycle (dummy pinch and estimation values)
orc.simple_solve(
    T_estimation_cs=10,
    T_estimation_hs=10,
    T_pinch_hs=5,
    T_pinch_cs=5
)

# Calculate thermal efficiency (ideal and real)
eta_ideal = orc.calc_thermal_efficiency(real=False)
eta_real = orc.calc_thermal_efficiency(real=True)
print("Thermal efficiency (ideal):", eta_ideal)
print("Thermal efficiency (real):", eta_real)

# Calculate mass flow from heat source
mass_flow_result = orc.calc_wf_mass_flow_from_heat_source(
    hs_mass_flow=89,
    p_hs=700*1000,
    p_cs=100000,
    hs_fluid_ref='Water',
    cs_fluid_ref='Air',
    real=True,
    print_table=True
)
wf_mass_flow  = mass_flow_result["wf"]["mass_flow"]
hs_mass_flow = mass_flow_result["hs"]["mass_flow"]
cs_mass_flow = mass_flow_result["cs"]["mass_flow"]
hs = {"h":mass_flow_result["hs"]["h"], "s":mass_flow_result["hs"]["s"]}
cs = {"h":mass_flow_result["cs"]["h"], "s":mass_flow_result["cs"]["s"]}

ex_dest = orc.calc_specific_exergy_destruction(
    T0=298.15,  # dead state temp [K]
    wf_mass_flow=wf_mass_flow,
    hs_mass_flow=hs_mass_flow,
    cs_mass_flow=cs_mass_flow,
    hs=hs,
    cs=cs,
    real=True
)

print("\nSpecific Exergy Destruction per Component:")
for comp, val in ex_dest.items():
    print(f"{comp:15s}: {val:.3f} W")

exer.plot_exergy_bar_chart(ex_dest)
exer.plot_exergy_pie_chart(ex_dest)

orc.print_state_table()
fig, ax = plt.subplots()
t_array_ideal,s_array_ideal = orc.get_Ts_arrays(real=False)
t_array_real,s_array_real = orc.get_Ts_arrays(real=True)
zoom1 = {
    "x0":   1.1155, "y0":24.26+273.15  ,
    "dx": 0.0091, "dy": 2.37,
    "zoom": 40, "loc": "lower left",
    "edgecolor": "green", "linestyle": "--", "linewidth": 1
}
zoom2 = {
    "x0": 1.7884, "y0": 105+273.15,
    "dx": 0.05, "dy": 7,
    "zoom": 10, "loc": "upper left",
    "edgecolor": "green", "linestyle": "--", "linewidth": 1
}
name_ideal = np.array(["0","1","2","23","3","4"])
name_real = np.array(["0_real","1_real","2_real","23_real","3_real","4_real"])
ts_plot=[(t_array_real,s_array_real),(t_array_ideal,s_array_ideal)]
# TSdiagram('R245fa',process_cycles=ts_plot)
TSdiagram("R245fa",
          process_cycles=ts_plot,
          ax=ax,
          zoom_insets=[zoom1, zoom2],
          name_ideal=name_ideal,
          name_real=name_real,figsize=(15,8))
#plt.savefig("figure.png", dpi=600)   # default is 100
plt.show()
# print("Mass flow result:", mass_flow_result)

# # Prepare dummy enthalpy/entropy data for exergy destruction
# hs = {
#     "h": [420e3, 410e3, 405e3],
#     "s": [1.2e3, 1.15e3, 1.1e3]
# }
# cs = {
#     "h": [100e3, 95e3, 90e3],
#     "s": [0.4e3, 0.38e3, 0.36e3]
# }

# exergy_destruction = orc.calc_specific_exergy_destruction(
#     T0=298.15,
#     wf_mass_flow=0.1,
#     hs_mass_flow=1.0,
#     cs_mass_flow=1.0,
#     hs=hs,
#     cs=cs,
#     real=True
# )
# print("Exergy destruction:", exergy_destruction)
# if __name__ == "__main__":
#     main()
