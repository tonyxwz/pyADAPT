import numpy as np
from scipy.integrate import odeint

# define model
def model1(x, t, p):
    # x: state variables on page 184
    # t: time variable to differntiation
    # p: parameter
    x_FC = x[0]
    x_CE_cyt = x[1]
    x_CE_ER = x[2]
    x_TG_cyt = x[3]
    x_TG_ER = x[4]
    x_TG_VLDL = x[5]
    x_C_VLDL = x[6]
    x_C_HDL = x[7]
    x_FFA = x[8]

    V_plasma = 1 # volume of blood plasma

    F_FC_prod = p[0]                # f1
    F_FC_met  = p[1] * x_FC         # f2
    F_CEfor_cyt = p[2] * x_FC       # f3
    F_CEdef_cyt = p[3] * x_CE_cyt   # f4
    F_CEfor_ER = p[4] * x_FC        # f5
    F_CEdef_ER = p[5] * x_CE_ER     # f6
    F_TGdnl_cyt = p[6]              # f7
    F_TGmet_cyt = p[7] * x_TG_cyt   # f8
    F_TGfor_cyt = p[8] * x_TG_ER    # f9
    F_TGdnl_ER = p[9]               # f10
    F_TGfor_ER = p[10] * x_TG_cyt   # f11
    F_FFA_upt = p[11] * x_FFA       # f12
    F_FFA_prod = p[12]              # f13
    F_VLDL_TG = p[13] * x_TG_ER     # f14
    F_VLDL_CE = p[14] * x_CE_ER     # f15
    F_TGupt_hep = p[15] * x_TG_VLDL # f16
    F_CEupt_hep = p[15] * x_C_VLDL  # f17
    F_TGupt_per = p[16] * x_TG_VLDL # f18
    F_CEupt_per = p[16] * x_C_VLDL  # f19
    F_CEfor_HDL = p[19]             # f20
    F_CEupt_HDL = p[20] * x_C_HDL   # f21
    F_TGhyd_hep = p[17] * x_TG_VLDL # f22
    F_TGhyd_per = p[18] * x_TG_VLDL # f23
    F_apoB_prod = p[21]             # f24

    dxdt = np.zeros(x.shape)

    dxdt[0] = F_FC_prod + F_CEdef_cyt + F_CEdef_ER - F_FC_met - F_CEfor_cyt - F_CEfor_ER
    dxdt[1] = F_CEfor_cyt - F_CEdef_cyt + V_plasma * (F_CEupt_hep + F_CEupt_HDL)
    dxdt[2] = F_CEfor_ER - F_CEdef_ER - F_VLDL_CE
    dxdt[3] = F_TGdnl_cyt + F_TGfor_cyt - F_TGfor_ER - F_TGmet_cyt + V_plasma * (F_FFA_upt / 3 + F_TGupt_hep + F_TGhyd_hep)
    dxdt[4] = F_TGdnl_ER + F_TGfor_ER - F_TGfor_cyt - F_VLDL_TG
    dxdt[5] = F_VLDL_TG / V_plasma - F_TGupt_hep - F_TGupt_per - F_TGhyd_hep - F_TGhyd_per
    dxdt[6] = F_VLDL_CE / V_plasma - F_CEupt_hep - F_CEupt_per
    dxdt[7] = F_CEfor_HDL - F_CEupt_HDL
    dxdt[8] = F_FFA_prod - F_FFA_upt
    
    return dxdt
