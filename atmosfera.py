import numpy as np
from scipy.interpolate import interp1d
#questa funzione mi serve per estrapolare una funzione che mi calcoli il termine sigma dell'ozono per ogni lunghezza d'onda

class Atmosfera:
    def __init__(self,
     enable_rayleigh: bool=True,
     enable_ozone: bool=True, 
     n_aria: float = 1.00029, 
     N_molecole: float = 2.504e25, 
     s_zenith: float = 8000.0,
     R_T: float = 6.371e6,
     N_O3: float = 8.06e22,
     percorso_file_O3: str | None = None,
     colonna_temperatura_ozono: int = 3):
        self.enable_rayleigh = enable_rayleigh
        self.enable_ozone = enable_ozone
        self.n_aria = n_aria
        self.N_molecole = N_molecole
        self.s_zenith = s_zenith
        self.R_T = R_T
        self.N_O3 = N_O3
        self.percorso_file_O3 = percorso_file_O3
        self.colonna_temperatura_ozono = colonna_temperatura_ozono
        self._o3_interp = None
        if self.enable_ozone and self.percorso_file_O3 is not None:
            self._load_ozone_data(self.percorso_file_O3, self.colonna_temperatura_ozono)

    def _load_ozone_data(self, percorso_file: str, colonna_temp: int):
        dati = np.loadtxt(percorso_file, comments=('#',';','*'))
        wavelength_m = dati[:,0] * 1e-9
        sigma_m2 = dati[:,colonna_temp] * 1e-4
        self._o3_interp = interp1d(wavelength_m, sigma_m2, kind="linear", bounds_error=False, fill_value=0.0)

    def percorso_ottico(self, angolo_zenitale_deg: float):
        theta = np.deg2rad(angolo_zenitale_deg)
        cos_t = np.cos(theta)
        return (np.sqrt((self.R_T*cos_t)**2 +2.0 * self.R_T * self.s_zenith +self.s_zenith**2) -self.R_T * cos_t)

    def coefficiente_rayleigh(self, wavelength: np.ndarray):
        num = 8.0 * (np.pi ** 3) *((self.n_aria**2 - 1.0) ** 2)
        den = 3.0 * (wavelength ** 4) * self.N_molecole
        return num/den

    def spessore_ottico(self,  wavelength: np.ndarray, angolo_zenitale_deg: float):
        s = self.percorso_ottico(angolo_zenitale_deg)
        tau = np.zeros_like(wavelength, dtype=float)
        if self.enable_rayleigh:
            tau += self.coefficiente_rayleigh(wavelength) * s
        if self.enable_ozone and self._o3_interp is not None:
            fattore_massa_aria = s / self.s_zenith
            sigma_o3 = self._o3_interp(wavelength)
            tau += sigma_o3 * (self.N_O3 * fattore_massa_aria)
        return tau

    def N_obs(self, flusso_N0: np.ndarray, wavelength: np.ndarray, angolo_zenitale_deg: float):
        tau = self.spessore_ottico(wavelength, angolo_zenitale_deg)
        return flusso_N0 * np.exp(-tau)
    
