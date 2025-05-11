#### Multiphase Flow - Mechanistic Flow Pattern Prediction - Mauricio Salkeld 11824210

import math

# Given Data
q_oil = 15  # m³/day
DI = 8.9
d = DI# cm
A = math.pi * (DI / 100)**2 / 4
Psc = 1 # bar (Surface Pressure)
Tlength = 3000  # m (total tubing length)
q = 30  # m3/day (flow rate)
wc = 0.5  # Water cut
API = 25  # Oil API gravity
Pintake = 250  # bar (Intake pressure)
Pbp = 200  # bar (Bubble point pressure)
sggas = 0.8  # Specific gravity of gas
Tsc = 16  # °C (Surface temperature)
Tin = 106  # °C (Inlet temperature)
R = 8.3144 # Gas constant
GOR = 100
Rs_3000 = GOR

# Geothermal and Pressure Gradients
pressure_gradient = 0.083  # bar/m
geothermal_gradient = 30 / 1000  # °C/m

# Defining depth sections
depth_sections = []
start_depth = Tlength
step = 300  # Each section decreases by 300m

for i in range(10):
    end_depth = start_depth - step
    if end_depth < 0:
        end_depth = 0  # Ensure no negative depths
    depth_sections.append((start_depth, end_depth))
    start_depth = end_depth

# Print sections to verify
for i, (start, end) in enumerate(depth_sections, 1):
    print(f"Section {i}: {start}m to {end}m")

# Calculate Oil Density (γ_oil)
gamma_oil = 141.5 / (API + 131.5)
gamma = gamma_oil  # alias for compatibility
rho_oil = gamma_oil * 1000  # Convert to kg/m3

# Water Density (Assumed standard value)
rho_water = 1000  # kg/m3

# Liquid Mixture Density
rho_L = wc * rho_water + (1 - wc) * rho_oil

# Calculate Bubble Point Depth and Temperature
bubble_point_depth = (Pbp - Psc) / pressure_gradient  # m
T_bp = geothermal_gradient * bubble_point_depth + Tsc  # °C

# Calculate Formation Volume Factor at Bubble Point
Rs_bp = GOR  # Rs at bubble point
if gamma_oil < 0.876:
    C1, C2, C3 = 2.622e-3, 1.100e-5, 1.337e-9
else:
    C1, C2, C3 = 2.626e-3, 1.751e-5, -1.811e-8

Bo_bp = 1 + (C1 * Rs_bp) + ((254.7 * (T_bp + 273.15)) / gamma_oil - (236.7 * (T_bp + 273.15)) - (73580 / gamma_oil) + 68380) * (C2 + C3 * Rs_bp)

# Oil Density at Bubble Point
T_bp_K = T_bp + 273.15
rho_o_bp = (1000 * gamma_oil + 1.224 * sggas * Rs_bp) / Bo_bp

# Oil viscosity at or below bubble point (Beggs-Robinson, using Kelvin)
def mu_dead_oil(T_K, gamma_o):
    T_R = T_K * 9/5  # Convert K to °R
    x = (T_R - 460) ** (-1.163) * math.exp(13.108 - 6.591 / gamma_o)
    mu_od = (10**(x - 3) - 0.001) / 1000  # Convert from cP to Pa·s
    return mu_od

def mu_oil_sat(T_K, Rs):
    mu_od = mu_dead_oil(T_K, gamma_oil)
    A = 10.715 * (5.615 * Rs + 100)**(-0.515)
    B = 5.44 * (5.615 * Rs + 150)**(-0.338)
    mu_os = A * mu_od**B
    return mu_os

mu_bp = mu_oil_sat(T_bp_K, Rs_bp)

# Oil viscosity above bubble point
def mu_gas(P_psia, T_R, SGg, z=1.0):
    Mg = 28.967 * SGg
    X = 3.5 + 986 / T_R + 0.001 * Mg
    Y = 2.4 - 0.2 * X
    K = ((9.4 + 0.02 * Mg) * T_R**1.5) / (209 + 19 * Mg + T_R)
    rho_g = (1 / 62.428) * (28.967 * SGg * P_psia) / (z * 10.732 * T_R)  # in g/cm³
    mu_g = K * math.exp(X * rho_g**Y) / 10000  # Convert to Pa·s from cP
    return mu_g

def mu_oil_undersat(P_kPa, mu_ob, Pb_kPa):
    m = 0.263 * P_kPa**1.187 * math.exp(-11.513 - 1.302e-5 * P_kPa)
    return mu_ob * (P_kPa / Pb_kPa)**m

# Define oil density formulas upfront
def rho_o_formula_undersat(P, Rs, T_K):
    return rho_o_bp * math.exp(((28.1 * Rs) + (30.6 * T_K) - (1180 * sggas) + (1784 / gamma_oil) - 10910) / (P * 1e5) * (P - Pbp))

def rho_o_formula_sat(Rs, Bo):
    return (1000 * gamma_oil + 1.224 * sggas * Rs) / Bo

# Saturated Rs calculation based on pressure and temperature
def Rs_saturated(P_kPa, T_K):
    yg = 1.225 + 0.00164 * T_K - 1.769 / gamma_oil
    return sggas * (P_kPa / (519.7 * (10**yg)))**1.204

# Function to calculate section properties
def calculate_section(depth, Rs=None):
    P = (Psc + pressure_gradient * depth) * 1e5  # Convert bar to Pa
    P_kPa = P / 1000
    T = Tsc + geothermal_gradient * depth  # °C
    T_K = T + 273.15  # Convert to Kelvin

    if P <= Pbp * 1e5:
        Rs = Rs_saturated(P_kPa, T_K)
    elif Rs is None:
        Rs = GOR

    Bo = 1 + (C1 * Rs) + ((254.7 * T_K) / gamma_oil - (236.7 * T_K) - (73580 / gamma_oil) + 68380) * (C2 + C3 * Rs)

    if P > Pbp * 1e5:
        rho_o = rho_o_formula_undersat(P, Rs, T_K)
        mu_o = mu_oil_undersat(P / 1000, mu_bp, Pbp * 100)
    else:
        rho_o = rho_o_formula_sat(Rs, Bo)
        mu_o = mu_oil_sat(T_K, Rs)

    Mr_gas = (sggas * 28.97) / 1000
    rho_gas = (P * Mr_gas) / (R * T_K)

    T_R = T_K * 9/5
    P_psia = P / 6894.76
    mu_g = None
    if P <= Pbp * 1e5:
        mu_g = mu_gas(P_psia, T_R, sggas)

    return P / 1e5, T, Bo, rho_o, Mr_gas, rho_gas, mu_o, mu_g

# Calculate for key depths
P_3000, T_3000, Bo_3000, rho_o_3000, _, _, mu_o_3000, mu_g_3000 = calculate_section(3000, Rs_3000)
P_2700, T_2700, Bo_2700, rho_o_2700, _, _, mu_o_2700, mu_g_2700 = calculate_section(2700, Rs_3000)
P_2400, T_2400, Bo_2400, rho_o_2400, _, _, mu_o_2400, mu_g_2400 = calculate_section(2400, Rs_3000)
P_2100, T_2100, Bo_2100, rho_o_2100, Mr_gas_2100, rho_gas_2100, mu_o_2100, mu_g_2100 = calculate_section(2100)

# Section 4: Calculate for 1800 m
P_1800, T_1800, Bo_1800, rho_o_1800, Mr_gas_1800, rho_gas_1800, mu_o_1800, mu_g_1800 = calculate_section(1800)

# Rs at 1800 m
Rs_1800 = Rs_saturated(P_1800 * 100, T_1800 + 273.15)

# Flow calculations at 1800 m
dq_gas_1800 = (Rs_bp - Rs_1800) * q_oil  # m³/day
q_gas_1800 = dq_gas_1800
q_gas_1800_m3s = q_gas_1800 / (24 * 3600)  # m³/s
v_sg_1800 = q_gas_1800_m3s / A

# Superficial liquid velocity at 1800 m
v_sl_1800 = (q_oil / (24 * 3600)) / A

# Terminal velocity at 1800 m
sigma_l = 0.03
sin_phi = 1
g = 9.81
v_t_1800 = 3.1 * ((sigma_l * g * sin_phi * (rho_o_1800 - rho_gas_1800)) / (rho_gas_1800 ** 2)) ** 0.25

# Rs at 2100 m
Rs_2100 = Rs_saturated(P_2100 * 100, T_2100 + 273.15)

# Flow calculations

dq_gas_2100 = (Rs_bp - Rs_2100) * q_oil  # m³/day
q_gas_2100 = dq_gas_2100
q_gas_2100_m3s = q_gas_2100 / (24 * 3600)

v_sg_2100 = q_gas_2100_m3s / A

# # Terminal velocity at 2100m
# sigma_l = 0.03
# sin_phi = 1
# g = 9.81
# v_t_2100 = 3.1 * ((sigma_l * g * sin_phi * (rho_o_2100 - rho_gas_2100)) / (rho_gas_2100 ** 2)) ** 0.25

# # Entrainment fraction f_E calculation
# term1 = (1e4 * v_sg_2100 * mu_g_2100) / sigma_l
# term2 = (rho_gas_2100 / rho_L) ** 0.5
# exponent = -0.125 * ((term1 * term2) - 1.5)
# exponent2 = math.exp(exponent)
# #f_E_2100 = (1 - exponent2)

# # Superficial liquid velocity [m/s]
# v_sl_2100 = (q_oil / (24 * 3600)) / A

# Liquid core fraction
# lambda_Lc_2100 = (v_sl_2100 * f_E_2100) / (v_sg_2100 + v_sl_2100 * f_E_2100)
# N_Re_f = (rho_L * v_sl_2100 * d * (1 - f_E_2100)) / mu_o_2100
# mu_o_2100f_SLf_f = 64 / N_Re_f
# N_Re_SL = (rho_L * v_sl_2100 * d) / mu_o_2100 
# f_SL = 64 / N_Re_SL

# # Modified composite density (Lockhart-Martinelli)
# rho_c_2100 = rho_gas_2100 * (1 - lambda_Lc_2100) + rho_L * lambda_Lc_2100
# mu_c_2100 = mu_g_2100 * (1 - lambda_Lc_2100) + mu_o_2100 * lambda_Lc_2100
# v_sc_2100 = v_sg_2100 + v_sl_2100 * f_E_2100
# N_Re_sc_2100 = (rho_c_2100 * v_sc_2100 * d) / mu_c_2100
# f_sc_2100 = 0.184 * N_Re_sc_2100**(-0.2)
# dp_dL_sc_2100 = - (f_sc_2100 * rho_c_2100 * v_sc_2100**2) / (2 * d)
# dp_dL_SL_2100 = - (f_SL * rho_L * v_sl_2100**2) / (2 * d)
# X_M_sq_2100 = dp_dL_SL_2100 / dp_dL_sc_2100
# Y_M_2100 = ((rho_L - rho_c_2100) * g) / abs(dp_dL_sc_2100)




# Print summary
print("--- Section 1: 2700 to 3000 m ---")
print(f"3000m - P: {P_3000:.2f} bar, T: {T_3000:.2f} °C, Bo: {Bo_3000:.4f}, rho_o: {rho_o_3000:.2f}, mu_o: {mu_o_3000:.6f}, mu_g: {mu_g_3000}")
print(f"2700m - P: {P_2700:.2f} bar, T: {T_2700:.2f} °C, Bo: {Bo_2700:.4f}, rho_o: {rho_o_2700:.2f}, mu_o: {mu_o_2700:.6f}, mu_g: {mu_g_2700}")
print("--- Section 2: 2400 to 2700 m ---")
print(f"2400m - P: {P_2400:.2f} bar, T: {T_2400:.2f} °C, Bo: {Bo_2400:.4f}, rho_o: {rho_o_2400:.2f}, mu_o: {mu_o_2400:.6f}, mu_g: {mu_g_2400}")
print("--- Section 3: 2100 to 2400 m ---")
print(f"2100m - P: {P_2100:.2f} bar, T: {T_2100:.2f} °C, Bo: {Bo_2100:.4f}, rho_o: {rho_o_2100:.2f}, mu_o: {mu_o_2100:.6f}, Rs: {Rs_2100:.2f}, Mr_gas: {Mr_gas_2100:.4f}, rho_gas: {rho_gas_2100:.2f}, mu_g: {mu_g_2100:.6f}")
print(f"q_oil: {q_oil:.2f} m³/day, q_gas_2100: {q_gas_2100:.2f} m³/day, v_sg_2100: {v_sg_2100:.4f} m/s")

# Section 4 summary
print("--- Section 4: 1800 to 2100 m ---")
print(f"1800m - P: {P_1800:.2f} bar, T: {T_1800:.2f} °C, Bo: {Bo_1800:.4f}, rho_o: {rho_o_1800:.2f}, mu_o: {mu_o_1800:.6f}, mu_g: {mu_g_1800}, v_t: {v_t_1800:.4f} m/s, q_oil: {q_oil:.2f} m³/day, q_gas_1800: {q_gas_1800:.2f} m³/day, v_sg_1800: {v_sg_1800:.4f} m/s")

# Section 5: 1500 to 1800 m
P_1500, T_1500, Bo_1500, rho_o_1500, Mr_gas_1500, rho_gas_1500, mu_o_1500, mu_g_1500 = calculate_section(1500)
Rs_1500 = Rs_saturated(P_1500 * 100, T_1500 + 273.15)
dq_gas_1500 = (Rs_bp - Rs_1500) * q_oil
q_gas_1500 = dq_gas_1500
q_gas_1500_m3s = q_gas_1500 / (24 * 3600)
v_sg_1500 = q_gas_1500_m3s / A
v_sl_1500 = (q_oil / (24 * 3600)) / A
v_t_1500 = 3.1 * ((sigma_l * g * sin_phi * (rho_o_1500 - rho_gas_1500)) / (rho_gas_1500 ** 2)) ** 0.25
print("--- Section 5: 1500 to 1800 m ---")
print(f"1500m - P: {P_1500:.2f} bar, T: {T_1500:.2f} °C, Bo: {Bo_1500:.4f}, rho_o: {rho_o_1500:.2f}, mu_o: {mu_o_1500:.6f}, mu_g: {mu_g_1500}, v_t: {v_t_1500:.4f} m/s, q_oil: {q_oil:.2f} m³/day, q_gas_1500: {q_gas_1500:.2f} m³/day, v_sg_1500: {v_sg_1500:.4f} m/s")

# Section 6: 1200 to 1500 m
P_1200, T_1200, Bo_1200, rho_o_1200, Mr_gas_1200, rho_gas_1200, mu_o_1200, mu_g_1200 = calculate_section(1200)
Rs_1200 = Rs_saturated(P_1200 * 100, T_1200 + 273.15)
dq_gas_1200 = (Rs_bp - Rs_1200) * q_oil
q_gas_1200 = dq_gas_1200
q_gas_1200_m3s = q_gas_1200 / (24 * 3600)
v_sg_1200 = q_gas_1200_m3s / A
v_sl_1200 = (q_oil / (24 * 3600)) / A
v_t_1200 = 3.1 * ((sigma_l * g * sin_phi * (rho_o_1200 - rho_gas_1200)) / (rho_gas_1200 ** 2)) ** 0.25
print("--- Section 6: 1200 to 1500 m ---")
print(f"1200m - P: {P_1200:.2f} bar, T: {T_1200:.2f} °C, Bo: {Bo_1200:.4f}, rho_o: {rho_o_1200:.2f}, mu_o: {mu_o_1200:.6f}, mu_g: {mu_g_1200}, v_t: {v_t_1200:.4f} m/s, q_oil: {q_oil:.2f} m³/day, q_gas_1200: {q_gas_1200:.2f} m³/day, v_sg_1200: {v_sg_1200:.4f} m/s")

# Section 7: 900 to 1200 m
P_900, T_900, Bo_900, rho_o_900, Mr_gas_900, rho_gas_900, mu_o_900, mu_g_900 = calculate_section(900)
Rs_900 = Rs_saturated(P_900 * 100, T_900 + 273.15)
dq_gas_900 = (Rs_bp - Rs_900) * q_oil
q_gas_900 = dq_gas_900
q_gas_900_m3s = q_gas_900 / (24 * 3600)
v_sg_900 = q_gas_900_m3s / A
v_sl_900 = (q_oil / (24 * 3600)) / A
v_t_900 = 3.1 * ((sigma_l * g * sin_phi * (rho_o_900 - rho_gas_900)) / (rho_gas_900 ** 2)) ** 0.25
print("--- Section 7: 900 to 1200 m ---")
print(f"900m - P: {P_900:.2f} bar, T: {T_900:.2f} °C, Bo: {Bo_900:.4f}, rho_o: {rho_o_900:.2f}, mu_o: {mu_o_900:.6f}, mu_g: {mu_g_900}, v_t: {v_t_900:.4f} m/s, q_oil: {q_oil:.2f} m³/day, q_gas_900: {q_gas_900:.2f} m³/day, v_sg_900: {v_sg_900:.4f} m/s")

# Section 8: 600 to 900 m
P_600, T_600, Bo_600, rho_o_600, Mr_gas_600, rho_gas_600, mu_o_600, mu_g_600 = calculate_section(600)
Rs_600 = Rs_saturated(P_600 * 100, T_600 + 273.15)
dq_gas_600 = (Rs_bp - Rs_600) * q_oil
q_gas_600 = dq_gas_600
q_gas_600_m3s = q_gas_600 / (24 * 3600)
v_sg_600 = q_gas_600_m3s / A
v_sl_600 = (q_oil / (24 * 3600)) / A
v_t_600 = 3.1 * ((sigma_l * g * sin_phi * (rho_o_600 - rho_gas_600)) / (rho_gas_600 ** 2)) ** 0.25
print("--- Section 8: 600 to 900 m ---")
print(f"600m - P: {P_600:.2f} bar, T: {T_600:.2f} °C, Bo: {Bo_600:.4f}, rho_o: {rho_o_600:.2f}, mu_o: {mu_o_600:.6f}, mu_g: {mu_g_600}, v_t: {v_t_600:.4f} m/s, q_oil: {q_oil:.2f} m³/day, q_gas_600: {q_gas_600:.2f} m³/day, v_sg_600: {v_sg_600:.4f} m/s")

# Section 9: 300 to 600 m
P_300, T_300, Bo_300, rho_o_300, Mr_gas_300, rho_gas_300, mu_o_300, mu_g_300 = calculate_section(300)
Rs_300 = Rs_saturated(P_300 * 100, T_300 + 273.15)
dq_gas_300 = (Rs_bp - Rs_300) * q_oil
q_gas_300 = dq_gas_300
q_gas_300_m3s = q_gas_300 / (24 * 3600)
v_sg_300 = q_gas_300_m3s / A
v_sl_300 = (q_oil / (24 * 3600)) / A
v_t_300 = 3.1 * ((sigma_l * g * sin_phi * (rho_o_300 - rho_gas_300)) / (rho_gas_300 ** 2)) ** 0.25
print("--- Section 9: 300 to 600 m ---")
print(f"300m - P: {P_300:.2f} bar, T: {T_300:.2f} °C, Bo: {Bo_300:.4f}, rho_o: {rho_o_300:.2f}, mu_o: {mu_o_300:.6f}, mu_g: {mu_g_300}, v_t: {v_t_300:.4f} m/s, q_oil: {q_oil:.2f} m³/day, q_gas_300: {q_gas_300:.2f} m³/day, v_sg_300: {v_sg_300:.4f} m/s")

# Section 10: 0 to 300 m (surface)
# Override PVT at surface: dead oil at Tsc, Psc
P_0 = Psc
T_0 = Tsc
Bo_0 = 1.0
rho_o_0 = rho_oil
mu_o_0 = mu_dead_oil(T_0 + 273.15, gamma_oil)
mu_g_0 = None
dq_gas_0 = 0.0
q_gas_0 = 0.0
v_sg_0 = 0.0
v_sl_0 = (q_oil / (24*3600)) / A
v_t_0 = 0.0
print("--- Section 10: 0 to 300 m ---")
print(f"0m  - P: {P_0:.2f} bar, T: {T_0:.2f} °C, Bo: {Bo_0:.4f}, rho_o: {rho_o_0:.2f}, mu_o: {mu_o_0:.6f}, mu_g: {mu_g_0}, v_t: {v_t_0:.4f} m/s, q_oil: {q_oil:.2f} m³/day, q_gas_0: {q_gas_0:.2f} m³/day, v_sg_0: {v_sg_0:.4f} m/s")
