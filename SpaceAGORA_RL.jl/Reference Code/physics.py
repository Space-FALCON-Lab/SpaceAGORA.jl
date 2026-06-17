import math
import numpy as np
import scipy
from scipy import special
def am(self, Area_SA, Area_SC, S):
    alpha = math.pi / 2
    sigma = 0.9

    def normalcoefficient(S, aoa, sigma):
        CN = 1 / (S ** 2) * ((((2 - sigma) / math.pi ** 0.5) * S * math.sin(aoa) + sigma / 2.0) * math.exp(
            -(S * math.sin(aoa)) ** 2) + (
                                     (2.0 - sigma) * (
                                     (S * math.sin(aoa)) ** 2.0 + 0.5) + sigma / 2.0 * math.pi ** 0.5 * (
                                             S * math.sin(aoa))) * (1.0 + math.erf(S * math.sin(aoa))))
        return CN

    def axialcoefficient(S, aoa, sigma):
        CA = ((sigma * math.cos(aoa)) / (math.pi ** 0.5 * S)) * (
                math.exp(-(S * math.sin(aoa)) ** 2) + math.pi ** 0.5 * (S * math.sin(aoa)) * (
                1.0 + math.erf(S * math.sin(aoa))))
        return CA

    ## Solar Panels:
    CN_sa = normalcoefficient(S, alpha, sigma)
    CA_sa = axialcoefficient(S, alpha, sigma)
    CL_sa = CN_sa * math.cos(alpha) - CA_sa * math.sin(alpha)
    CD_sa = CA_sa * math.cos(alpha) + CN_sa * math.sin(alpha)

    ## Spacecraft
    CN_sc = normalcoefficient(S, math.pi * 0.5, sigma)
    CA_sc = axialcoefficient(S, math.pi * 0.5, sigma)
    CL_sc = CN_sc * math.cos(math.pi * 0.5) - CA_sc * math.sin(math.pi * 0.5)
    CD_sc = CA_sc * math.cos(math.pi * 0.5) + CN_sc * math.sin(math.pi * 0.5)

    CD_body = (CD_sa * Area_SA + CD_sc * Area_SC) / (Area_SA + Area_SC)
    CL_body = (CL_sa * Area_SA + CL_sc * Area_SC) / (Area_SA + Area_SC)

    return CL_body, CD_body

def ht(self,v,T,rho):
    R = self.p.R
    aoa = math.pi/2
    thermal_accomodation_factor = 1

    T_p = T
    T_w = T
    S = (v / (2 * R*T) ** 0.5)
    gamma = self.p.gamma

    first_term = np.multiply(1e-4 * thermal_accomodation_factor * R * T_p * (np.sqrt(R * T_p / (2 * math.pi))),
                             rho)
    term_a = np.exp(-np.square(np.multiply(S,np.sin(aoa)))) + np.multiply(np.sqrt(math.pi) * (np.multiply(S,np.sin(aoa))),
                            (1 + scipy.special.erf(np.multiply(S,np.sin(aoa)))))
    term_b = np.multiply((S ** 2 + (gamma) / (gamma - 1) - (gamma + 1) / (2 * (gamma - 1)) * (T_w / T_p)),term_a)

    heat_rate = np.multiply(term_b- 0.5 * np.exp(
                        -(np.square(np.multiply(S,np.sin(aoa))))),first_term )     # W/cm^2

    return heat_rate