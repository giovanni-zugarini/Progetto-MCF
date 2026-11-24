import numpy as np
from scipy.constants import h, c, k

class CorpoNero:

    def __init__(self, nome: str, Temperatura: float, raggio: float, distanza: float):
        self.nome = nome
        self.Temperatura = Temperatura
        self.raggio = raggio
        self.distanza = distanza

    def angolo_solido(self):
        return np.pi * (self.raggio/self.distanza) ** 2

    def radianza_spettrale(self, wavelength: np.ndarray):
        esponente = (h*c) / (wavelength * k * self.Temperatura)
        return (2.0 * h * c**2 /wavelength**5)/np.expm1(esponente)

    def densità_fotoni(self, wavelength: np.ndarray):
        E_gamma = (h*c)/wavelength
        return self.radianza_spettrale(wavelength)/E_gamma

    def toa_flusso_fotoni(self, wavelength: np.ndarray):
        return self.densità_fotoni(wavelength) * self.angolo_solido()
#l'ultima funzione mi serve per moltiplicare D (densità di fotoni) è una grandezza specifica per angolo solido, infatti nella 
# sua unità di misura è presente sr^-1, moltiplicando per l'angolo solido elimino la dipendenza dagli steradianti, ottenendo 
# il flusso incidente sull'atmosfera 