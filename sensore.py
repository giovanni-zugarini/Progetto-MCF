import numpy as np


class sensore:
    def __init__(self, nome: str, area: float, sensibilità: float, t_esposizione: float, dlambda: float):
       
        self.nome = nome
        self.area = area
        self.sensibilità = sensibilità
        self.t_esposizione = t_esposizione
        self.d_lamba = self.d_lamba

    def sensibilità(self, wavelength:np.ndarray):
        if callable(self.sensibilità):
            return self.sensibilità(wavelength)
        return np.full_like(wavelength, float(self.sensibilità, dtype=float))

    def conteggio_medio_atteso(self, n_obs_flux: np.ndarray, wavelength: np.ndarray)
        eps=self.sensibilità(wavelength)
        return n_obs_flux*self.dlambda*self.area*self.t_esposizione*eps
    
    def misura(self, n_obs_flux: np.ndarray, wavelenght: np.ndarray):
            N = self.conteggio_medio_atteso(n_obs_flux, wavelenght)
            K = self.rng.poisson(N)
            sigma = np.sqrt(N)
            SNR=np.divide(N, sigma, out=np.zeros_like(N), where=sigma>0)
            return {'atteso':N, 'misurato':K, 'incertezza':sigma, 'SNR':SNR}