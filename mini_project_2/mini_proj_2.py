# Mini Project 2

import numpy as np
import matplotlib.pyplot as plt
import photovoltaic as pv


###################################################################################################################
# diffusion length function
# input: doping level (Na or Nd), C_auger, teff, D
# output: diffusion length in cm
###################################################################################################################
def diffusion_length(doping, auger_coefficient, undoped_lifetime, diffusivity):
    # simplified formula for diffusion length from doping
    auger_lifetime = 1 / (auger_coefficient * doping ** 2)  # auger lifetime (s)
    lifetime = 1 / (1 / auger_lifetime + 1 / undoped_lifetime)  # lifetime = 1/SRH + 1/auger
    # print(lifetime)
    return np.sqrt(lifetime * diffusivity)  # diffusion length (cm)


###################################################################################################################
# photon flux at wavelength function
# input: spectral irradiance, wavelength
# output: photon flux
###################################################################################################################
def ph_flux(spec_irr, w_length):
    p_flux = spec_irr * (wavelength / (1240 * q_j))
    return p_flux


###################################################################################################################
# emitter resistance function
# input: finger spacing (cm), emitter sheet resistivity (ohm/ sqr)
# output: emitter resistance (ohm cm^2)
###################################################################################################################
def emitter_resistance(Sf, Rsheet):
    return Rsheet * Sf ** 2 / 12


###################################################################################################################
# sheet resistivity function
# input: doping, thickness
# output: sheet resistivity
###################################################################################################################
def sheet_resistivity(doping, thickness):
    mobility = pv.si.mob_masetti_phos(doping)
    # mobility = 300  # assume constant mobility
    return 1 / (q_j * doping * mobility * thickness)


#######################################################################################################################
# Part 1
###############################
# parameters for calculations #
###############################
# all for Si
ni = 8.6e9  # intrinsic carrier concentration cm-3
De = 15  # diffusivity emitter cm2/s
Db = 30  # diffusivity base cm2/s
Ne = 1e18  # doping emitter cm-3
Nb = 5e16  # doping base cm-3
Se = 500  # surface recombination emitter cm/s
Sb = 1000  # surface recombination base cm/s
We = 0.0001  # width emitter cm
Wb = 0.03  # width base cm
q_j = 1.6e-19  # charge of electron in joules
q_eV = 1  # charge of electron in eV
kT = 0.0259  # Boltzmann constant
# effective minority carrier lifetime calculations
t0 = 1e-3  # minority carrier lifetime sec
C_auger = 1.5e-30
inv_teff = C_auger * ni ** 2 + (1 / t0)
teff = (1 / inv_teff)  # effective minority carrier lifetime s/cm6
Le = diffusion_length(Ne, C_auger, teff, De)  # diffusion length in the emitter(cm)
Lb = diffusion_length(Nb, C_auger, teff, Db)  # diffusion length in the base(cm)
# Le = np.sqrt(De * teff)  # diffusion length emitter cm
# Lb = np.sqrt(Db * teff)  # diffusion length base cm
mu_e = (De * q_j) / kT  # mobility emitter side
mu_b = (Db * q_j) / kT  # mobility base side
data = np.loadtxt('10_nm_AM15G_absorption.txt')  # reads absorption coefficient and AM1.5G spectra from file
wavelength = data[:, 0]  # wavelength (nm) column numbering starts at 0
am15g = data[:, 1]  # AM15G spectrum (W/m2)
absorption_coefficient = data[:, 2]  # absorption coefficient (/cm)

# output 1 plot minority carrier lifetime as a function of doping
Ne_plot = np.linspace(ni, Ne, num=100)  # doping array
Le_plot = diffusion_length(Ne_plot, C_auger, teff, De)  # diffusion length array
t_plot = (Le_plot / De)  # minority carrier lifetime array
plt.figure('minority carrier lifetime vs doping')
plt.title('Minority Carrier Lifetime in Si')
plt.xlabel('Doping (cm^-3)')
plt.ylabel('Minority Carrier Lifetime (sec)')
plt.xscale('log')
plt.yscale('log')
plt.plot(Ne_plot, t_plot)
plt.show()

# output 2 plot absorption coefficient vs wavelength
plt.figure('absorption coefficient vs wavelength')
plt.title('Absorption Coefficient as a Function of Wavelength')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Absorption Coefficient (cm^-1)')
plt.yscale('log')
plt.plot(wavelength, absorption_coefficient)
plt.show()

# output 3 print values in readable format
print('All values for SI')
print('Intrinsic Carrier Concentration (cm^-3):', ni)
print('Minority Carrier Lifetime (sec):', t0)
print('Diffusivity in Emitter (cm^2/s):', De)
print('Diffusivity in Base (cm^2/s):', Db)
print('Doping in Emitter (cm^-3):', Ne)
print('Doping in Base (cm^-3):', Nb)
print('Surface Recombination in Emitter (cm/s:', Se)
print('Surface Recombination in Base (cm/s:', Sb)
print('Thickness of Emitter (cm):', We)
print('Thickness of Base (cm):', Wb)
print('Diffusion Length of Emitter (cm):', Le)
print('Diffusion Length of Base (cm):', Lb)
print('Mobility in Emitter (cm^2/V*s):', mu_e)
print('Mobility in Base (cm^2/V*s):', mu_b)
print('Effective Minority Carrier Lifetime (s/cm^6):', teff)

#######################################################################################################################
# Part 2
# question 1 calculate J0
J0e = pv.cell.J0_layer(We, Ne, De, Le, Se, ni)  # dark saturation current in the emitter (A/cm2)
J0b = pv.cell.J0_layer(Wb, Nb, Db, Lb, Sb, ni)  # dark saturation current in the base (A/cm2)
J0 = J0e + J0b
print('J0 emitter (A/cm²):', J0e)
print('J0 base (A/cm²):', J0b)
print('J0 total (A/cm²):', J0)

# question 2
# a
photon_flux = ph_flux(am15g, wavelength)
plt.figure('photon flux vs wavelength')
plt.title('Photon Flux at Each Wavelength')
plt.xlabel('wavelength (nm)')
plt.ylabel('Photon Flux (m^-2*s^-1')
plt.plot(wavelength, photon_flux)
plt.show()

# b
# quantum efficiency
IQEemitter = pv.cell.iqe_emitter(absorption_coefficient, We, Le, De, Se)
IQEbase = pv.cell.iqe_base(absorption_coefficient, We, Wb, Lb, Db, Sb)
IQEtotal = IQEemitter + IQEbase
plt.figure('IQE')  # start a new figure
plt.title('Quantum Efficiency as a Function of Wavelength')
plt.ylabel('Quantum Efficiency')
plt.xlabel('Wavelength (nm)')
plt.plot(wavelength, IQEemitter, label='IQE Emitter')
plt.plot(wavelength, IQEbase, label='IQE Base')
plt.plot(wavelength, IQEtotal, label='IQE Total')
plt.legend()
plt.show()

# c
# JL

JL_base = photon_flux * q_j * IQEbase  # A/m^2
JL_base = JL_base * 1e-4  # A/cm^2
JL_emitter = photon_flux * q_j * IQEemitter  # A/m^2
JL_emitter = JL_emitter * 1e-4  # A/cm^2
JL = photon_flux * q_j * IQEtotal  # A/m^2
JL = JL * 1e-4  # A/cm^2
JL_total = 0  # placeholder
JL_base_total = 0  # placeholder
JL_emitter_total = 0  # placeholder
# loops to sum over all wavelengths for JL totals
for i in JL:
    JL_total += i
for i in JL_base:
    JL_base_total += i
for i in JL_emitter:
    JL_emitter_total += i

plt.figure('JLi')  # start a new figure
plt.title('Total JL (A/cm²): 0.03894')
plt.ylabel('JL (A/cm²)')
plt.xlabel('Wavelength (nm)')
plt.plot(wavelength, JL)
plt.show()

print('light generated current density in emitter (A/cm²):', JL_emitter_total)
print('light generated current density in base (A/cm²):', JL_base_total)
print('total light generated current density (A/cm²):', JL_total)

# question 3
# calculate emitter series resistance
Sf = 0.2  # metal finger spacing at the front (cm)
Rsheet = sheet_resistivity(Ne, We)
Rseries = emitter_resistance(Sf, Rsheet)  # (ohm cm2) only consider the series resistance from the emitter
print('series resistance (ohm.cm²)', Rseries)
# plot the IV Curve
# At Voc generation balances with recombination
Voc = pv.cell.Voc(JL_total, J0)
# print('Voc', Voc)
# sweep from zero volts to Voc
voltage = np.linspace(0, Voc, 100)
current = pv.cell.I_cell(voltage, JL_total, J0)
Voc, Jsc, FF, Vmp, Jmp = pv.cell.cell_params(voltage, current)
# print(Voc, Isc, FF, Vmp, Jmp, Vmp*Jmp)
voltage_r = voltage - current * Rseries
Voc_r, Jsc_r, FF_r, Vmp_r, Jmp_r = pv.cell.cell_params(voltage_r, current)
# print(Voc_r, Jsc_r, FF_r, Vmp_r, Jmp_r, Pmp_r)
plt.figure('iv curve')
plt.title('IV Curve')
plt.ylabel('Current density (A/cm²)')
plt.xlabel('Voltage (V)')
# plt.plot(voltage, current)
plt.plot(voltage_r, current)
plt.show()

#######################################################################################################################
# Part 3
# all of the following values calculated in part 2 and printed here
# question 1
print('Jsc:', Jsc)
# question 2
Voc_calc = (0.0259/q_eV) * np.log(Jsc/J0)  # calculated instead of using param function
print('Voc:', Voc_calc)
# question 3
print('FF with series resistance', FF_r)
# question 4
print('Efficiency (%)', Vmp_r * Jmp_r * 1000)  # 1000 = 100/0.1 W/cm²

#######################################################################################################################
# Part 4
# all done including resistance
# optimize for efficiency by looping through all emitter width values at each doping values and running calcs
Ne_opt = np.linspace(1e10, 1e20, num=400)  # doping array
Le_opt = diffusion_length(Ne_opt, C_auger, teff, De)  # diffusion length array
Le_max = np.amax(Le_opt)  # max diffusion length
We_opt = np.linspace(0.00000001, Le_max, num=400)  # emitter width array
eff_opt = []  # placeholder for efficiency values
We_maxeff = 0  # placeholder
Ne_maxeff = 0  # placeholder
Voc_maxeff = 0  # placeholder
Jsc_maxeff = 0  # placeholder
FF_maxeff = 0  # placeholder
for i in Ne_opt:
    for j in We_opt:
        # find saturation current
        J0e_opt = pv.cell.J0_layer(j, i, De, Le, Se, ni)  # dark saturation current in the emitter (A/cm2)
        J0_opt = J0e_opt + J0b  # total dark saturation current
        # find quantum efficiency
        IQEemitter_opt = pv.cell.iqe_emitter(absorption_coefficient, j, Le, De, Se)
        IQEbase_opt = pv.cell.iqe_base(absorption_coefficient, j, Wb, Lb, Db, Sb)
        IQEtotal_opt = IQEemitter_opt + IQEbase_opt
        # find light generated current
        JL_opt = photon_flux * q_j * IQEtotal_opt  # A/m^2
        JL_opt = JL_opt * 1e-4  # A/cm^2
        JL_total_opt = 0  # placeholder
        # loops to sum over all wavelengths for JL total
        for k in JL_opt:
            JL_total_opt += k
        Rsheet_opt = sheet_resistivity(i, j)
        Rseries_opt = emitter_resistance(Sf, Rsheet_opt)  # (ohm cm2) only consider the resistance from the emitter
        Voc_opt = pv.cell.Voc(JL_total_opt, J0_opt)
        current_opt = pv.cell.I_cell(voltage, JL_total_opt, J0_opt)
        voltage_r_opt = voltage - current_opt * Rseries_opt
        Voc_r_opt, Jsc_r_opt, FF_r_opt, Vmp_r_opt, Jmp_r_opt = pv.cell.cell_params(voltage_r_opt, current_opt)
        eff = Vmp_r_opt * Jmp_r_opt * 1000  # efficiency
        if eff == 22.15074000449771:  # this is computed max eff in loop (line 271), change this if max eff changes
            Ne_maxeff = i  # getting Ne value at max eff
            We_maxeff = j  # getting We value at max eff
            Voc_maxeff = Voc_r_opt  # Voc at max eff
            Jsc_maxeff = Jsc_r_opt  # Jsc at max eff
            FF_maxeff = FF_r_opt  # FF at max eff
        eff_opt = np.append(eff_opt, eff)  # append to eff array for extraction
max_eff = np.amax(eff_opt)  # maximum efficiency value
print('Emitter doping value at maximum efficiency is (cm^-3):', Ne_maxeff)
print('Emitter width value at maximum efficiency is (cm):', We_maxeff)
print('Jsc value at maximum efficiency is (A/cm^2):', Jsc_maxeff)
print('Voc value at maximum efficiency is (V):', Voc_maxeff)
print('Fill Factor value at maximum efficiency is:', FF_maxeff)
print('Maximum efficiency with resistance is (%):', max_eff)
#######################################################################################################################
# Part 5
# analyze the losses in the solar cell
# change material parameters and plot against efficiency to find the largest effect
# doping vs efficiency
eff_adj_Ne_array = []  # placeholder
for i in Ne_opt:
    # find saturation current
    J0e_adj_Ne = pv.cell.J0_layer(We, i, De, Le, Se, ni)  # dark saturation current in the emitter (A/cm2)
    J0_adj_Ne = J0e_adj_Ne + J0b  # total dark saturation current
    # find quantum efficiency
    IQEemitter_adj_Ne = pv.cell.iqe_emitter(absorption_coefficient, We, Le, De, Se)
    IQEbase_adj_Ne = pv.cell.iqe_base(absorption_coefficient, We, Wb, Lb, Db, Sb)
    IQEtotal_adj_Ne = IQEemitter_adj_Ne + IQEbase_adj_Ne
    # find light generated current
    JL_adj_Ne = photon_flux * q_j * IQEtotal_adj_Ne  # A/m^2
    JL_adj_Ne = JL_adj_Ne * 1e-4  # A/cm^2
    JL_total_adj_Ne = 0  # placeholder
    # loops to sum over all wavelengths for JL total
    for k in JL_adj_Ne:
        JL_total_adj_Ne += k
    Rsheet_adj_Ne = sheet_resistivity(i, We)
    Rseries_adj_Ne = emitter_resistance(Sf, Rsheet_adj_Ne)  # (ohm cm2) only consider the resistance from the emitter
    Voc_adj_Ne = pv.cell.Voc(JL_total_adj_Ne, J0_adj_Ne)
    current_adj_Ne = pv.cell.I_cell(voltage, JL_total_adj_Ne, J0_adj_Ne)
    voltage_r_adj_Ne = voltage - current_adj_Ne * Rseries_adj_Ne
    Voc_r_adj_Ne, Jsc_r_adj_Ne, FF_r_adj_Ne, Vmp_r_adj_Ne, Jmp_r_adj_Ne = pv.cell.cell_params(voltage_r_adj_Ne, current_adj_Ne)
    eff_adj_Ne = Vmp_r_adj_Ne * Jmp_r_adj_Ne * 1000  # efficiency
    if eff_adj_Ne <= 0:
        eff_adj_Ne_array = np.append(eff_adj_Ne_array, 0)  # append to eff array
    else:
        eff_adj_Ne_array = np.append(eff_adj_Ne_array, eff_adj_Ne)  # append to eff array
plt.figure('doping vs eff')
plt.title('Doping vs Efficiency')
plt.ylabel('Efficiency (%)')
plt.xlabel('Doping Level (cm^-3)')
plt.xscale('log')
plt.plot(Ne_opt, eff_adj_Ne_array)
plt.show()
# emitter width vs efficiency
eff_adj_We_array = []  # placeholder
for i in We_opt:
    # find saturation current
    J0e_adj_We = pv.cell.J0_layer(i, Ne, De, Le, Se, ni)  # dark saturation current in the emitter (A/cm2)
    J0_adj_We = J0e_adj_We + J0b  # total dark saturation current
    # find quantum efficiency
    IQEemitter_adj_We = pv.cell.iqe_emitter(absorption_coefficient, i, Le, De, Se)
    IQEbase_adj_We = pv.cell.iqe_base(absorption_coefficient, i, Wb, Lb, Db, Sb)
    IQEtotal_adj_We = IQEemitter_adj_We + IQEbase_adj_We
    # find light generated current
    JL_adj_We = photon_flux * q_j * IQEtotal_adj_We  # A/m^2
    JL_adj_We = JL_adj_We * 1e-4  # A/cm^2
    JL_total_adj_We = 0  # placeholder
    # loops to sum over all wavelengths for JL total
    for k in JL_adj_We:
        JL_total_adj_We += k
    Rsheet_adj_We = sheet_resistivity(Ne, i)
    Rseries_adj_We = emitter_resistance(Sf, Rsheet_adj_We)  # (ohm cm2) only consider the resistance from the emitter
    Voc_adj_We = pv.cell.Voc(JL_total_adj_We, J0_adj_We)
    current_adj_We = pv.cell.I_cell(voltage, JL_total_adj_We, J0_adj_We)
    voltage_r_adj_We = voltage - current_adj_We * Rseries_adj_We
    Voc_r_adj_We, Jsc_r_adj_We, FF_r_adj_We, Vmp_r_adj_We, Jmp_r_adj_We = pv.cell.cell_params(voltage_r_adj_We, current_adj_We)
    eff_adj_We = Vmp_r_adj_We * Jmp_r_adj_We * 1000  # efficiency
    eff_adj_We_array = np.append(eff_adj_We_array, eff_adj_We)  # append to eff array for extraction
plt.figure('we vs eff')
plt.title('Emitter Width vs Efficiency')
plt.ylabel('Efficiency (%)')
plt.xlabel('Emitter Width (cm)')
plt.xscale('log')
plt.plot(We_opt, eff_adj_We_array)
plt.show()



