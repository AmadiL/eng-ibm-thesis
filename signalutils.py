def normalized(func):
    def wrapper(*args, **kwargs):
        return normalize(func(*args, **kwargs))
    return wrapper


def normalize(signal):
    return signal / max(abs(signal))


def find_zeros(signal):
    indices = []
    for i in range(len(signal)):
        if not i == 0 and signal[i] * signal[i - 1] < 0:
            zero_idx = i if abs(signal[i]) < abs(signal[i - 1]) else (i - 1)
            indices.append(zero_idx)
    return indices


def samples_to_ms(samples, fs):
    return 1e3 * samples / fs


def ms_to_samples(milliseconds, fs):
    return fs * milliseconds / 1e3
