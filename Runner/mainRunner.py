# This file is part of LAO-STO.
#
# Copyright (C) 2025 Julian Czarnecki
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# If you use this code for scientific research, please cite:
# J. Czarnecki et. al.,
# "Superconducting gap symmetry of 2DEG at (111)-oriented LaAlO3/SrTiO3 interface",
# arXiv:2508.05075 (2025).
# https://arxiv.org/abs/2508.05075

from RunnerClass import *
import re
import numpy as np

def runTemperatureDependence():
    runner = Runner()
    pathToT0 = os.path.join(SCRATCH_PATH, "STO-SC", "LAO-STO-E_Fermi_J_SC")

    n_sublattices = 2
    Sublat_param = ("discretization", "SUBLATTICES", n_sublattices)
    const_v_layer = 0.0e3
    V_layer_param = (
        "physical_params",
        "V_layer",
        [const_v_layer for _ in range(n_sublattices)],
    )
    n_subbands = 1
    Subband_param = ("discretization", "SUBBANDS", n_subbands)
    Subband_energies_param = ("physical_params", "Subband_energies", [0.0])
    # Fermi energy
    nml_name = "physical_params"
    param_name = "E_Fermi"
    Ef_min = -1.05e3
    Ef_max = -0.9e3
    Ef_steps = 75
    dE = abs(Ef_max - Ef_min) / Ef_steps
    Fermi_table = [(nml_name, param_name, Ef_min + i * dE) for i in range(Ef_steps + 1)]

    # J_SC
    # nml_name = "physical_params"
    # param_name = "J_SC"
    # J_sc = 50.0
    # J_sc_nml = (nml_name, param_name, J_sc)
    # J_min = 0.36e3
    # J_max = 0.37e3
    # J_steps = 1
    # dJ = abs(J_max - J_min) / J_steps
    # J_table = [(nml_name, param_name, J_min + i * dJ) for i in range(J_steps + 1)]

    # J_SC_NNN
    nml_name = "physical_params"
    param_name = "J_SC_NNN"
    J_sc_nnn = 75.0
    J_sc_nnn_nml = (nml_name, param_name, J_sc_nnn)


    nml_name = "physical_params"
    param_name = "T"
    T_min = 0.65
    T_max = 1
    T_steps = 7
    dT = abs(T_max - T_min) / T_steps
    T_table = [(nml_name, param_name, T_min + i * dT) for i in range(T_steps + 1)]

    read_gamma = ("self_consistency", "read_gamma_from_file", True)
    read_charge = ("self_consistency", "read_charge_from_file", True)

    for T in T_table:
        for Ef in Fermi_table:
            t0DirGamma = os.path.join(pathToT0, f"RUN_E_Fermi_{Ef[2]}_J_SC_NNN_{J_sc_nnn_nml[2]}", "OutputData", "Gamma_SC_final.dat")
            t0DirCharge = os.path.join(pathToT0, f"RUN_E_Fermi_{Ef[2]}_J_SC_NNN_{J_sc_nnn_nml[2]}", "OutputData", "Charge_dens_final.dat")
            if not os.path.exists(t0DirGamma):
                t0DirGamma = os.path.join(pathToT0, f"RUN_E_Fermi_{Ef[2]}_J_SC_NNN_{J_sc_nnn_nml[2]}", "OutputData", "Gamma_SC_iter.dat")
                t0DirCharge = os.path.join(pathToT0, f"RUN_E_Fermi_{Ef[2]}_J_SC_NNN_{J_sc_nnn_nml[2]}", "OutputData", "Charge_dens_iter.dat")

            gammaDirSetting = ("self_consistency", "path_to_gamma_start", t0DirGamma)
            chargeDirSetting = ("self_consistency", "path_to_charge_start", t0DirCharge)

            runner.runSlurmParamValue(
                paramValuePairs=[
                    T,
                    Ef,
                    J_sc_nnn_nml,
                    read_gamma,
                    read_charge,
                    gammaDirSetting,
                    chargeDirSetting
                ],
                runsDir="STO-SC/LAO-STO-E_Fermi_J_SC_NNN_T",
                material="STO",
                isAres=True,
            )

def configureAndRunPostprocessing():
    runner = Runner()
    pathToRuns = os.path.join(SCRATCH_PATH, "KTO-SC", "KTO-E_Fermi_J_SC_NNN")
    # Adding '' at the end to terminate path with /
    # directories = [
    #     os.path.join(pathToRuns, dir, "")
    #     for dir in os.listdir(pathToRuns)
    #     if re.match("RUN.*", dir)
    # ]
    directories = [
        os.path.join(pathToRuns, f"RUN_E_Fermi_{ef}.0_J_SC_NNN_350.0", "")
        for ef in [-80, -72, -60, -52, -40, -32, -20, -12, 0]]

    enable_sc = ("sc_gap_calculation", "enable_sc_gap_calc", True)
    enable_chern = ("chern_number_calculation", "enable_chern_number_calc", False)
    enable_dos = ("dos_calculation", "enable_dos_calc", True)
    enable_gamma_k = ("gamma_k_calculation", "enable_gamma_k_calc", False)
    for dir in directories:
        nmlDirectorySC = ("sc_gap_calculation", "path_to_run_dir_sc_gap", dir)
        nmlDirectoryChern = (
            "chern_number_calculation",
            "path_to_run_dir_chern_number",
            dir,
        )
        nmlDirectoryDos = ("dos_calculation", "path_to_run_dir_dos", dir)
        nmlDirectoryGammaK = ("gamma_k_calculation", "path_to_run_dir_gamma_k", dir)
        runner.runSlurmPostprocessing(
            dir,
            [enable_sc, nmlDirectorySC, enable_dos, nmlDirectoryDos],
            machine="default"
        )


def configureAndRunSc():
    runner = Runner()

    #Fermi energy
    nml_name = "physical_params"
    param_name = "E_Fermi"
    Ef_min = 0.00e3
    Ef_max = 0.15e3
    Ef_steps = 50
    dE = abs(Ef_max - Ef_min) / Ef_steps
    Fermi_table = [(nml_name, param_name, Ef_min + i * dE) for i in range(Ef_steps + 1)]
    #Fermi_table = [(nml_name, param_name, -0.06e3)]

    # J_SC
    nml_name = "physical_params"
    param_name = "J_SC_tensor"
    J_SC_val = [0.0, 0.0, 0.0, 0.0, 0.0, 170, 0.0, 0.0, 0.0, 0.0, 170.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    J_SC = (nml_name, param_name, J_SC_val)

    #for phi in phi_table:
    for Ef in Fermi_table:
        runner.runSlurmParamValue(
            paramValuePairs=[
                Ef,
                J_SC,
                #U_nml,
                #V_nml
                #J_NNN
                #V_layer_param,
                #Sublat_param,
                #Subband_param,
                #Subband_energies_param,
            ],
            runsDir="STO-SC/STO-E_Fermi_J_SC",
            material="STO",
            machine="helios",
        )


def runDosFitting():
    runner = Runner()
    gammaMax = 0.8
    dGamma = 0.05

    nml_name = "physical_params"
    param_name = "E_Fermi"
    Ef_min = 0.02e3
    Ef_max = 0.06e3
    Ef_steps = 2
    dE = abs(Ef_max - Ef_min) / Ef_steps
    Fermi_table = [(nml_name, param_name, Ef_min + i * dE) for i in range(Ef_steps + 1)]

    nml_name = "physical_params"
    param_name = "B_field"

    Bz = 30 #mT
    B_field = (nml_name, param_name, [0, 0, Bz / 1e3])


    for Ef in Fermi_table:
        dosFitDict = {
            "runsDir": os.path.join(SCRATCH_PATH, 'KTO-SC', 'NicolasDosFitting', f'J1_{Bz}mT_Ef_{Ef[2]}'),
            "dosExpPath": os.path.join(HOME_PATH, 'NicolasFit', f'J1_exp_{Bz}mT.dat'),
            "eMax": 0.3,
            #"bounds": [[3.034e-01 - dGamma, 3.034e-01 + dGamma], [4.384e-01 + dGamma, 4.384e-01 - dGamma], [-4.408e-01 - dGamma, -4.408e-01 + dGamma]]
            "bounds": [[-gammaMax, gammaMax], [-gammaMax, gammaMax], [-gammaMax, gammaMax]]
        }

        runner.runSlurmDosFitter(fitConfig=dosFitDict,
                                    paramValuePairs=[Ef, B_field],
                                    material="KTO",
                                    isAres=True)


def runMockedOutputPostprocessing():
    runner = Runner()

    nml_name = "physical_params"
    param_name = "E_Fermi"
    Ef_min = 0.05e3
    Ef_max = 0.1e3
    Ef_steps = 2
    dE = abs(Ef_max - Ef_min) / Ef_steps
    Fermi_table = [(nml_name, param_name, Ef_min + i * dE) for i in range(Ef_steps + 1)]

    symmetriesWeightsDict = {"nearest": {
                                r"$A_1^{(1)}$": 0.1,
                                r"$A_1^{(2)}$": 0.4,
                                },
                            "next": {
                                r"$E_2^{(1)}$": 0.5,
                                r"$E_2^{(2)}$": 0.5,
                                },
                            }
    gammaAmplitudesDict = {"nearest": np.complex128(0.),
                            "next": np.complex128(0.1)}

    for Ef in Fermi_table:
        runner.runSlurmMockedOutputPostprocessing(
            paramValuePairs=[Ef],
            paramValuePairsPost=[],
            gammaAmplitudesDict=gammaAmplitudesDict,
            symmetriesWeightsDict=symmetriesWeightsDict,
            runsDir="STO-SC/LAO-STO-E1/",
            material="STO",
            machine="default",
        )

def getJTensor():
    runner = Runner()
    jCouplings = {("S", "S"): 350.0 * 2,}
    nmlInput = runner.createJTensorTable(jCouplings)
    print(nmlInput)

def main():
    #runTemperatureDependence()
    configureAndRunSc()
    #configureAndRunPostprocessing()
    #runDosFitting()
    #runMockedOutputPostprocessing()
    #getJTensor()

if __name__ == "__main__":
    main()
