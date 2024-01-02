#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 17:34:12 2022
@author: ben
"""

# from sklearn.metrics import r2_score
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from pyCATHY.plotters import cathy_plots as pltCT


def fit_VGP(sai, teta, SWRfunc, solver="lm", bvals=[(0, 0, 1, 0), (1, 1, 100, 100)]):
    """Fit VG parameters using Least square"""

    # Solver= 'dogbox' # select 'lm', 'trf' or 'dogbox'
    if solver == "lm":
        init_vals = [
            0.05,
            0.45,
            1,
            1,
        ]  # Four parameters tetar,tetas,n or lamda, alpha or beta
        fitP, covar = curve_fit(SWRfunc, sai, teta, p0=init_vals)  # Call solver
    else:
        init_vals = [
            0.05,
            0.45,
            1,
            1,
        ]  # Four parameters--tetar,tetas,n or lamda, alpha or sib
        fitP, covar = curve_fit(
            SWRfunc, sai, teta, p0=init_vals, method=solver, bounds=bvals
        )  # Solver

    tetamodel = SWRfunc(sai, fitP[0], fitP[1], fitP[2], fitP[3])  # fitted SWR function
    error = [np.absolute(covar[i][i]) ** 0.5 for i in range(len(fitP))]
    R2 = r2_score(teta, tetamodel)

    return fitP, error, R2


def fit_Archie(df):
    popt, pcov = curve_fit(
        Archie_sat2rho,
        xdata=df["soil"].to_numpy(),
        ydata=df["ert"].to_numpy(),
        p0=[1.0, 1.0],
    )
    return (popt, pcov)


def fit_Archie_TsoetAl(df):
    """To obtain best-fit estimates of Archie parameters, a straight
    line is fitted for log 10 (S) and log 10 (ρ S ) using the least-squares criterion. The fitting rou-
    tine returns the covariance structure of the model estimates, which can be used to de-
    termine the 68% confidence interval (1 standard deviation) of the model estimates.
    """
    pass


def PH2SW_VG_model(
    psi, tetar, tetas, n, psib, ax=None
):  # Model for function for VG model
    """
    Translate Pressure heads to Saturation Water using VG model
    """
    # theta(psi) relationship
    # https://www.sciencedirect.com/science/article/pii/S0098300421001898?via%3Dihub#appsec1
    alpha = (
        1 / psib
    )  # Paramètre de Van Genuchten (Unité de pression −1, : m^−1 ou Pa^−1 )
    m = 1 - (1 / n)
    teta = tetar + ((tetas - tetar) * (1 / (1 + (alpha * psi) ** n)) ** m)
    if ax is not None:
        theta_VG_df = pd.DataFrame(np.array([psi, teta]).T)
        theta_VG_df.columns = ["psi", "theta"]
        pltCT.plot_VGP(theta_VG_df)
        return teta, ax

    return teta


def PH2SW_BC_model(psi, tetar, tetas, lamda, sib):  # Model for function BC model
    """
    Translate Pressure heads to Saturation Water using BC model
    """
    # https://www.sciencedirect.com/science/article/pii/S0098300421001898?via%3Dihub#appsec1
    teta = np.zeros(len(psi))
    for i in range(len(psi)):
        if psi[i] > sib:
            teta[i] = tetar + ((tetas - tetar) * (psi[i] / sib) ** -lamda)
        else:
            teta[i] = tetas
    return teta


def SW2PH_VG_model(SW):
    """
    Translate Saturation Water to Pressure heads using VG model
    """
    PH = []  # pressure heads
    return PH


def Archie_sat2rho(eff_sat, rho_sat, n):
    print(rho_sat, n)
    rho = rho_sat * eff_sat ** (-n)
    # plt.plot(eff_sat, rho, 'o')
    # plt.pause(1)
    return rho


def Archie_rho2sat(rho, rFluid, porosity, a=1.0, m=2.0, n=2.0):
    """
    rho: resistivity
    𝑆𝑤 : water saturation
    𝜙: the porosity of the soil
    𝜎_{𝑤} is the conductivity of the pore fluid
    𝑎, 𝑚, and 𝑛 are empirically derived parameters
    𝑎 is the tortuosity factor
    𝑚 is the cementation exponent
    𝑛 is the saturation exponent
    Returns
    -------
    𝑆𝑤 : water saturation
    """
    return (rho / rFluid / porosity ** (-m)) ** ((-1 / n))


#%% ---------------------------------------------------------------------------
#############  predict soil physical properties from literature ###############
## ---------------------------------------------------------------------------


def predict_unsat_soil_hydr_param(data=[[20, 20, 60]]):
    """
    https://github.com/usda-ars-ussl/rosetta-soil

    The Rosetta pedotransfer function predicts five parameters for
    the van Genuchten model of unsaturated soil hydraulic properties

    - theta_r : residual volumetric water content
    - theta_s : saturated volumetric water content
    - log10(alpha) : retention shape parameter [log10(1/cm)]
    - log10(n) : retention shape parameter
    - log10(ksat) : saturated hydraulic conductivity [log10(cm/d)]

    Model Code

    - 2	sa, si, cl (SSC)
    - 3	SSC, bulk density (BD)
    - 4	SSC, BD, th33
    - 5	SSC, BD, th33, th1500

    With

    - sa, si, cl are percentages of sand, silt and clay
    - BD is soil bulk density (g/cm3)
    - th33 is the soil volumetric water content at 33 kPa
    - th1500 is the soil volumetric water content at 1500 kPa

    Parameters
    ----------
    data : TYPE, optional
        Sand, silt, and clay %. The default is [[20,20,60]].

    Returns
    -------
    VGP_fit : TYPE
        DESCRIPTION.

    """
    from rosetta import rosetta, SoilData

    # data = [[20,20,60,1.45]] # Sand, silt, and clay
    soildata = SoilData.from_array(data)
    mean, stdev, codes = rosetta(3, soildata)  # rosetta version 3

    VGP_predict_names = [
        "theta_r",
        "theta_s",
        "log10(alpha)",
        "log10(n)",
        "log10(ksat)",
    ]
    VGP_predict = {}
    for i, VGP in enumerate(VGP_predict_names):
        # print(i,VGP,mean[:, i])
        VGP_predict[VGP] = abs(float(mean[:, i]))
        VGP_predict[VGP + "_stdev"] = float(stdev[:, i])

        # array  |
        # column | parameter
        # -----------------
        #    0   | theta_r, residual water content (cm3/cm3)
        #    1   | theta_s, saturated water content (cm3/cm3)
        #    2   | log10(alpha), van Genuchten 'alpha' parameter (1/cm)
        #    3   | log10(npar), van Genuchten 'n' parameter
        #    4   | log10(Ksat), saturated hydraulic conductivity (cm/day)

    VGP_predict_CATHY = {}
    VGP_predict_CATHY["POROS"] = 1- VGP_predict["theta_s"]
    VGP_predict_CATHY["VGNCELL"] =  10 ** VGP_predict["log10(n)"]
    VGP_predict_CATHY["VGRMCCELL"] = VGP_predict["theta_r"]
    # VGP_predict_CATHY["VGPSATCELL"] = VGP_predict["theta_s"] #(1 / (10 ** VGP_predict["log10(alpha)"])) # * 1e-3  #  1/alpha
    VGP_predict_CATHY["VGPSATCELL"] = (1 / (10 ** VGP_predict["log10(alpha)"])) # * 1e-3  #  1/alpha
    VGP_predict_CATHY["PERMX"] = 10 ** (VGP_predict["log10(ksat)"]) * (1e-2 / 86400)
    VGP_predict_CATHY["PERMY"] = 10 ** (VGP_predict["log10(ksat)"]) * (1e-2 / 86400)
    VGP_predict_CATHY["PERMZ"] = 10 ** (VGP_predict["log10(ksat)"]) * (1e-2 / 86400)
    VGP_predict_CATHY["ELSTOR"] = 1e-5

    # import numpy as np

    # np.exp(1.169) * (1e-3 / 86400)
    # np.exp(2.808) * (1e-3 / 86400)

    return VGP_predict_CATHY


def past_authors(selec=0):
    """data = [weil,busato,botto_clay,botto_sand]"""

    # 4.05278E-05    4.05278E-05    4.05278E-05   0.001   0.41  PERMX PERMY PERMZ ELSTOR POROS st1
    # 2.28   0.057  -0.0806           VGN,VGRMC,VGPSAT
    # 4.05278E-05    4.05278E-05    4.05278E-05   0.001   0.41  2.28   0.057  -0.0806 PERMX PERMY PERMZ ELSTOR POROS
    # PERMX PERMY  PERMZ  ELSTOR POROS,VGNCELL,VGRMCCELL,VGPSATCELL

    # weil et al
    # ---------------------------------------------------------
    weil = [1.88e-04, 1.88e-04, 1.88e-04, 1.00e-05, 0.55, 1.46, 0.15, 0.03125]
    busato = [4.05278e-05, 4.05278e-05, 4.05278e-05, 0.001, 0.41, 2.28, 0.057, -0.0806]

    botto_sand = [1e-4, 1e-4, 1e-4, 5e-3, 0.58, 0.065, 0.070, 1.88]
    botto_clay = [1e-7, 1e-7, 1e-7, 5e-3, 0.40, 0.067, 0.017, 1.40]

    data = [weil, busato, botto_clay, botto_sand]
    VGP_predict_CATHY = {}
    VGP_predict_CATHY["PERMX"] = data[selec][0]
    VGP_predict_CATHY["PERMY"] = data[selec][0]
    VGP_predict_CATHY["PERMZ"] = data[selec][0]
    VGP_predict_CATHY["ELSTOR"] = data[selec][3]
    VGP_predict_CATHY["POROS"] = data[selec][4]
    VGP_predict_CATHY["VGNCELL"] = data[selec][5]
    VGP_predict_CATHY["VGRMCCELL"] = data[selec][6]
    VGP_predict_CATHY["VGPSATCELL"] = data[selec][7]


    return VGP_predict_CATHY


def Twarakavi_etal_2009(soilTexture):
    """
    An objective analysis of the dynamic nature of field capacity
    https://doi.org/10.1029/2009WR007944
    """

    df_Twarakavi_etal_2009_VGN = []

    return df_Twarakavi_etal_2009_VGN


def Leij_etal_1996():
    """

    @book{leij1996unsoda,
      title={The UNSODA unsaturated soil hydraulic database: user's manual},
      author={Leij, Feike J},
      volume={96},
      number={95},
      year={1996},
      publisher={National Risk Management Research Laboratory, Office of Research and~…}
    }

    The UNSODA unsaturated soil hydraulic database:
    user’s manual, volume 96. National Risk Management Research Laboratory,
    office of Research and Development, US Environmental Protection Agency.


    """

    # Sand 0.045 0.43 0.145 2.68 712.80
    # Loamy sand 0.057 0.41 0.124 2.28 350.16
    # Sandy loam 0.065 0.41 0.075 1.89 106.08
    # Loam 0.078 0.43 0.036 1.56 24.96
    # Silt 0.034 0.46 0.016 1.37 6.00
    # Silt loam 0.067 0.45 0.020 1.41 10.80
    # Sandy clay loam 0.100 0.39 0.059 1.48 31.44
    # Clay loam 0.095 0.41 0.019 1.31 6.24
    # Silty clay loam 0.089 0.43 0.010 1.23 1.68
    # Sandy clay 0.100 0.38 0.027 1.23 2.88
    # Silty clay 0.070 0.36 0.005 1.09 0.48
    # Clay 0.068 0.38 0.008 1.09 4.80

    return


def VGN_Nielsen_1985():
    """

    https://www.ars.usda.gov/ARSUserFiles/20360500/pdf_pubs/P0871.pdf

    On Describing and Predicting the Hydraulic Properties of Unsaturated Soils
    """

    pass


def Carsel_Parrish_1988(soilTexture=None):
    """
    Transformation of VGN parameters to make them normaly distributed

    https://doi.org/10.1029/WR024i005p00755

    S, sand; SL, sandy loam; LS, loamy sand; SIL, silt loam; SI, silt; C, clay

    Returns
    -------
    None.

    """

    # Ks hydraulic Conductivity  Ks (cm/h)
    # theta_r, residual water content ()
    # theta_s, saturated water content ()
    # alpha, retention parameter (1/cm)
    # N, water retention model parameter
    # ------------------------------------
    hydraulic_variable = ["theta_s", "theta_r", "Ks", "alpha", "N"]

    # statistical description
    # -----------------------------------
    description = ["xmean", "sd", "CV", "n"]

    # soil texture 1
    # ---------------------------------
    # 'xmean', 'sd','CV','n'
    # theta_s      ..     ..    ..  ..
    # theta_r      ..     ..    ..  ..
    # Ks           ..     ..    ..  ..
    # alpha        ..     ..    ..  ..
    # N            ..     ..    ..  ..

    # soil texture 2
    # ---------------------------------
    # 'xmean', 'sd','CV','n'
    # theta_s      ..     ..    ..  ..
    # theta_r      ..     ..    ..  ..
    # Ks           ..     ..    ..  ..
    # alpha        ..     ..    ..  ..
    # N            ..     ..    ..  ..

    soil_texture = ["C", "S", "SI", "SIL", "CL"]

    if soilTexture is not None:
        if soilTexture not in soil_texture:
            print("soil texture non existing")

    table = [
        [
            [0.38, 0.09, 24.1, 400],
            [0.068, 0.034, 49.9, 353],
            [0.20, 0.42, 210.3, 114],
            [0.008, 0.012, 160.3, 400],
            [1.09, 0.09, 7.9, 400],
        ],
        [
            [0.43, 0.06, 15.1, 246],
            [0.045, 0.010, 22.3, 246],
            [29.70, 15.60, 52.4, 246],
            [0.145, 0.029, 20.3, 246],
            [2.68, 0.29, 20.3, 246],
        ],
        [
            [0.46, 0.11, 17.4, 82],
            [0.034, 0.010, 29.8, 82],
            [0.25, 0.33, 129.9, 88],
            [0.016, 0.007, 45.0, 82],
            [1.37, 0.05, 3.3, 82],
        ],
        [
            [0.45, 0.08, 18.7, 1093],
            [0.067, 0.015, 21.6, 1093],
            [0.45, 1.23, 275.1, 1093],
            [0.020, 0.012, 64.7, 1093],
            [1.41, 0.12, 8.5, 1093],
        ],
        [
            [0.41, 0.09, 22.4, 364],
            [0.095, 0.010, 10.1, 363],
            [0.26, 0.70, 267.2, 345],
            [0.019, 0.015, 77.9, 363],
            [1.31, 0.09, 7.2, 364],
        ],
    ]

    dict_Carsel_Parrish_1988_VGN = {}
    for i, st in enumerate(soil_texture):
        dict_Carsel_Parrish_1988_VGN[st] = {}
        for j, hv in enumerate(hydraulic_variable):
            # dict_Carsel_Parrish_1988_VGN[st][hv] = {}
            conversion = 1
            if "Ks" in hv:
                conversion = 1e-3 / 3600
            elif "alpha" in hv:
                conversion = 1 / 1e-3

            dict_Carsel_Parrish_1988_VGN[st][hv] = {
                "xmean": table[i][j][0] * conversion,
                "sd": table[i][j][1],
                "CV": table[i][j][2],
                "n": table[i][j][3],
            }

    # 29.7*1e-3

    df_Carsel_Parrish_1988_VGN = (
        pd.DataFrame.from_dict(dict_Carsel_Parrish_1988_VGN).stack().to_frame()
    )
    df_Carsel_Parrish_1988_VGN = pd.DataFrame(
        df_Carsel_Parrish_1988_VGN[0].values.T.tolist(),
        index=df_Carsel_Parrish_1988_VGN.index,
    )
    df_Carsel_Parrish_1988_VGN.index.names = ["hydraulic_variable", "soil_texture"]
    df_Carsel_Parrish_1988_VGN = df_Carsel_Parrish_1988_VGN.reorder_levels(
        ["soil_texture", "hydraulic_variable"]
    ).sort_index()

    # transformations
    # ------------------------------------------------------------------------
    hydraulic_variable = ["Ks", "theta_r", "alpha", "N"]
    soil_texture = ["C", "S", "SL", "SIL", "CL"]

    if soilTexture is not None:
        if soilTexture not in soil_texture:
            print("soil texture non existing")

    A = [
        [0, 0, 0, 0.9],
        [0, 0, 0, 1.5],
        [0, 0, 0, 1.35],
        [1, 0, 0, 1.35],
        [0, 0, 0, 1],
    ]

    B = [
        [5, 0.15, 0.15, 1.4],
        [70, 0.1, 0.25, 4],
        [30, 0.11, 0.25, 3],
        [15, 0.11, 0.15, 2],
        [7.5, 0.13, 0.15, 1.6],
    ]

    transf_type = [
        ["SB", "SUt", "SBt", "LNt"],
        ["SB", "LN", "SB", "LN"],
        ["SB", "SB", "LN", "SB"],
        ["LN", "SB", "LN", "SB"],
        ["SBtt", "SU", "LN", "SB"],
    ]

    dict_Carsel_Parrish_1988_transf = {}
    for i, st in enumerate(soil_texture):
        dict_Carsel_Parrish_1988_transf[st] = {}
        for j, hv in enumerate(hydraulic_variable):
            dict_Carsel_Parrish_1988_transf[st][hv] = {
                "A": A[i][j],
                "B": B[i][j],
                "transf_type": transf_type[i][j],
            }

    df_Carsel_Parrish_1988_transf = (
        pd.DataFrame.from_dict(dict_Carsel_Parrish_1988_transf).stack().to_frame()
    )
    df_Carsel_Parrish_1988_transf = pd.DataFrame(
        df_Carsel_Parrish_1988_transf[0].values.T.tolist(),
        index=df_Carsel_Parrish_1988_transf.index,
    )
    df_Carsel_Parrish_1988_transf.index.names = ["hydraulic_variable", "soil_texture"]
    df_Carsel_Parrish_1988_transf = df_Carsel_Parrish_1988_transf.reorder_levels(
        ["soil_texture", "hydraulic_variable"]
    ).sort_index()

    # between 0.36 (close random packing) and 0.40 (loose random packing) [Allen, 1985]

    if soilTexture is not None:
        df_Carsel_Parrish_1988_VGN = df_Carsel_Parrish_1988_VGN.xs(
            soilTexture, level="soil_texture"
        )
        df_Carsel_Parrish_1988_transf = df_Carsel_Parrish_1988_transf.xs(
            soilTexture, level="soil_texture"
        )

        VGP_predict_CATHY = {}

        # if soilTexture == 'S':
        #     # [Allen, 1985]
        #     VGP_predict_CATHY['POROS']=  0.36
        # if soilTexture == 'C':
        #     # [Allen, 1985]
        #     VGP_predict_CATHY['POROS']=  0.36

        VGP_predict_CATHY["POROS"] = df_Carsel_Parrish_1988_VGN["xmean"]["theta_s"]
        VGP_predict_CATHY["VGNCELL"] = df_Carsel_Parrish_1988_VGN["xmean"]["N"]
        VGP_predict_CATHY["VGRMCCELL"] = df_Carsel_Parrish_1988_VGN["xmean"]["theta_r"]
        VGP_predict_CATHY["VGPSATCELL"] = (
            1 / df_Carsel_Parrish_1988_VGN["xmean"]["alpha"]
        )
        VGP_predict_CATHY["PERMX"] = df_Carsel_Parrish_1988_VGN["xmean"]["Ks"]
        VGP_predict_CATHY["PERMY"] = df_Carsel_Parrish_1988_VGN["xmean"]["Ks"]
        VGP_predict_CATHY["PERMZ"] = df_Carsel_Parrish_1988_VGN["xmean"]["Ks"]

    return df_Carsel_Parrish_1988_VGN, df_Carsel_Parrish_1988_transf, VGP_predict_CATHY


def FAO_56_chapter7():

    pass


def Melo_Lier_2021():
    """
    https://doi.org/10.1016/j.jhydrol.2021.126952

    """

    feddes_CATHY = []

    return feddes_CATHY


def Feddes_litterature(crop=None):
    """
    Collect Feddes values from literature
    https://www.sciencedirect.com/science/article/pii/S1878029613002776?via%3Dihub

    Returns
    -------
    None.

    """

    dict_feddes = {}

    crops = ["Corn", "Soybeans"]

    #  anaerobiosis point (θan)
    # incipient water stress (θref);
    # he wilting point (θwp)

    # 'PCANA':0, 'PCREF':-4, 'PCWLT':-150,'ZROOT':0.2, 'PZ':1, 'OMGC':1,

    h1 = [-0.15, -0.10]  # s point
    h2 = [-0.30, -0.25]  # anaerobiosis point PCANA
    h3 = [-5.00, -8.00]  # incipient water stress PCREF
    h4 = [-80.00, -160.00]  # wilting point pressure head PCWLT

    # cols_name=['crop','h1','h2','h3','h4','ref']

    refs = [
        "https://doi.org/10.1016/0022-1694(83)90045-8",
        "https://doi.org/10.1016/0378-3774(94)90041-8",
    ]

    for i, c in enumerate(crops):
        dict_feddes[c] = {
            "h1": h1[i],
            "PCANA": h2[i],
            "PCREF": h3[i],
            "PCWLT": h4[i],
            "ref": refs[i],
        }
    df_feddes = pd.DataFrame(dict_feddes)

    if crop is not None:
        df_feddes = df_feddes[crop]

    return df_feddes


def Leij_etal_1996():
    """
    Description of soil physical parameters

    .. seealso::  Leij, F. J., W. J. Alves, and M. T. van Genuchten, 1996. The UNSODA unsaturated soil hydraulic
                    database: user’s manual, volume 96. National Risk Management Research Laboratory, Oﬃce of
                    Research and Development, US Environmental Protection Agency.

    Returns
    -------
    None.

    """
