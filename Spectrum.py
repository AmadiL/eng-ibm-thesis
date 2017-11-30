from scipy.fftpack import *

class Spectrum:
    """
    Spectrum of a signal.
    """
    def __init__(self, y, t, N, dt):
        self.y = rfft(y)
        self.f = rfftfreq(N, dt)
        self.N = len(self.f)
        self.fs = 1 / dt