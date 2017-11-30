import os
import wfdb
import numpy as np
from scipy.signal import *
from Spectrum import Spectrum
import matplotlib.pyplot as plt
import peakutils
from signalutils import *
import scipy.fftpack as fft
import bisect

SAMPLE_DATA_DIR = os.path.dirname(os.path.realpath(__file__)) + '/sampledata'


class ECGAnalyzer:
    """
    Class for ECG analysis.
    """
    def __init__(self):
        self.sample_data_path = SAMPLE_DATA_DIR + '/patient001/s0010_re'
        self.wfdb_record = wfdb.rdsamp(self.sample_data_path)

        self.fs = self.wfdb_record.fs # Hz
        self.dt = 1 / self.fs
        self.nyquist_f = self.fs/2 # Hz
        self.lowpass_fc = 15 # Hz
        self.highpass_fc = 2 # Hz

        self.lowpass_order = 1
        self.highpass_order = 2
        self.integration_window_length = int(0.15 * self.fs)

        self.ecg_signal = [np.array(i, dtype=np.float32) for i in zip(*self.wfdb_record.p_signals)][0]
        self.length = self.ecg_signal.__len__()
        self.t = np.arange(self.length) / self.fs

        self.no_drift_ecg_signal = None
        self.filtered_ecg_signal = None
        self.differentiated_ecg_signal = None
        self.squared_ecg_signal = None
        self.integrated_ecg_signal = None

    def calculate(self):
        # Bandpass filtration
        self.no_drift_ecg_signal = self._highpass(self.ecg_signal)
        self.filtered_ecg_signal = self._lowpass(self.no_drift_ecg_signal)
        # Differentiation
        self.differentiated_ecg_signal = self._derivative(self.filtered_ecg_signal)
        # Squaring
        self.squared_ecg_signal = self._squaring(self.differentiated_ecg_signal)
        # Moving-Window Integration
        self.integrated_ecg_signal = self._moving_window_integral(self.squared_ecg_signal)
        # plt.plot(self.t, self.differentiated_ecg_signal)
        # plt.plot(self.t, self.filtered_ecg_signal)
        # plt.show(block=True)
        self._detection_laguna_thakor(self.differentiated_ecg_signal)

    def _detection_pan_tompkins(self, signal):
        peaks = peakutils.indexes(signal, thres=0)
        last_peak_idx = 0
        running_estimate_signal_peak = 0
        running_estimate_noise_peak = 0
        first_peak_threshold = 0

        qrs_indices = np.empty(0, dtype=int)
        noise_indices = np.empty(0, dtype=int)

        for peak_idx, peak_val in zip(peaks, signal[peaks]):
            if peak_val > first_peak_threshold:
                qrs_indices = np.append(qrs_indices, peak_idx)
                running_estimate_signal_peak = 0.125 * peak_val + 0.875 * running_estimate_signal_peak
            else:
                noise_indices = np.append(noise_indices, peak_idx)
                running_estimate_noise_peak = 0.125 * peak_val + 0.875 * running_estimate_noise_peak
            first_peak_threshold = running_estimate_noise_peak \
                                   + 0.25 * (running_estimate_signal_peak - running_estimate_noise_peak)
        return qrs_indices, noise_indices

    def _detection_laguna_thakor(self, signal):
        peak_indecies = peakutils.indexes(np.abs(signal), thres=0)
        zero_indecies = find_zeros(signal)
        peak_n = max(abs(signal[:self.fs*2]))
        threshold_n = 0.8 * peak_n
        RR_interval_avg = None
        # RR_intervals = np.empty(0, dtype=int)

        # R_indecies = np.empty(0, dtype=int)
        # S_indecies = np.empty(0, dtype=int)
        # Q_indecies = np.empty(0, dtype=int)
        # QRS_start_indecies = np.empty(0, dtype=int)
        # P_indecies = np.empty(0, dtype=int)
        # T_indecies = np.empty(0, dtype=int)
        # T_end_indecies = np.empty(0, dtype=int)
        # QT_intervals = np.empty(0, dtype=int)
        # QT_peak_intervals = np.empty(0, dtype=int)
        # QT_corrected_intervals = np.empty(0, dtype=int)
        # QT_peak_corrected_intervals = np.empty(0, dtype=int)

        indecies = {
            'Q': np.empty(0, dtype=int),
            'R': np.empty(0, dtype=int),
            'S': np.empty(0, dtype=int),
            'QRS_start': np.empty(0, dtype=int),
            'P': np.empty(0, dtype=int),
            'T': np.empty(0, dtype=int),
            'T_end': np.empty(0, dtype=int),
        }

        intervals = {
            'RR': np.empty(0, dtype=int),
            'QT': np.empty(0, dtype=int),
            'QTP': np.empty(0, dtype=int),
            'QTc': np.empty(0, dtype=int),
            'QTPc': np.empty(0, dtype=int)
        }

        for (idx, peak_idx), peak_val in zip(enumerate(peak_indecies), signal[peak_indecies]):
            """
            Iterate over peak values of differentiated signal
            idx - index of array with peak indecies
            peak_idx - index of the peak -> peak_indecies[idx]
            peak_val - value of the differentiated signal at peak_idx -> signal[peak_idx]
            """
            if abs(peak_val) > threshold_n and (len(indecies['R']) == 0 or (peak_idx - indecies['R'][-1]) > 0.2 * self.fs):
                # If peak value is over threshold and one of conditions is met:
                # 1) None R peak was found yet
                # 2) Peak is more than 200ms away from last detected R peak

                # Find Q, R and S peaks + adjust threshold
                Q_idx, R_idx, S_idx, threshold_n = self._LT_QRS_detect(threshold_n=threshold_n,
                                                                       idx=idx,
                                                                       peak_idx=peak_idx,
                                                                       peak_val=peak_val,
                                                                       peak_indecies=peak_indecies,
                                                                       zero_indecies=zero_indecies,
                                                                       signal=signal)
                if len(indecies['R']):
                    # If at least one R peak was recorded calculated RR interval and its average
                    # Cannot be executed at first R peak
                    RR_interval_avg, RR_interval = self._LT_RR_interval_calc(indecies=indecies,
                                                                             intervals=intervals,
                                                                             RR_interval_avg=RR_interval_avg,
                                                                             R_idx=R_idx)
                # Record R and S peaks
                indecies['R'] = np.append(indecies['R'], R_idx)
                indecies['S'] = np.append(indecies['S'], S_idx)

                # Detect QRS beginning
                QRS_start_idx = self._LT_QRS_start_detect(R_idx=R_idx,
                                                          Q_idx=Q_idx,
                                                          indecies=indecies,
                                                          peak_indecies=peak_indecies,
                                                          signal=signal)
                # Record QRS beginning
                indecies['QRS_start'] = np.append(indecies['QRS_start'], QRS_start_idx)

                if len(indecies['T_end']):
                    # If at least one end of T wave was recorded detect P peak
                    P_idx = self._LT_P_detect(indecies=indecies, QRS_start_idx=QRS_start_idx, signal=signal)
                    # Record P peak
                    indecies['P'] = np.append(indecies['P'], P_idx)

                if RR_interval_avg:
                    # If RR interval average was calculated detect T peak and end of T wave
                    T_end_idx, T_idx = self._LT_T_detect(RR_interval_avg=RR_interval_avg,
                                                         R_idx=R_idx,
                                                         indecies=indecies,
                                                         zero_indecies=zero_indecies,
                                                         signal=signal)
                    # Record T peak and end of T wave
                    indecies['T'] = np.append(indecies['T'], T_idx)
                    indecies['T_end'] = np.append(indecies['T_end'], T_end_idx)
                    # Calculate QT, QTP, QTc and QTPc intervals [ms]
                    QT_intervals = self._LT_intervals_calc(T_end_idx=T_end_idx,
                                                           T_idx=T_idx,
                                                           QRS_start_idx=QRS_start_idx,
                                                           RR_interval=RR_interval)
                    QT_interval, QT_peak_interval, QT_corrected_interval, QT_peak_corrected_interval = QT_intervals
                    # Record QT, QTP, QTc and QTPc intervals [ms]
                    intervals['QT'] = np.append(intervals['QT'], QT_interval)
                    intervals['QTP'] = np.append(intervals['QTP'], QT_peak_interval)
                    intervals['QTc'] = np.append(intervals['QTc'], QT_corrected_interval)
                    intervals['QTPc'] = np.append(intervals['QTPc'], QT_peak_corrected_interval)

        plt.figure(0)
        plt.plot(self.t, self.filtered_ecg_signal)
        plt.axhline(0, color='black')
        plt.plot(self.t[indecies['Q']], self.filtered_ecg_signal[indecies['Q']], 'ro')
        plt.plot(self.t[indecies['R']], self.filtered_ecg_signal[indecies['R']], 'go')
        plt.plot(self.t[indecies['S']], self.filtered_ecg_signal[indecies['S']], 'yo')
        plt.plot(self.t[indecies['QRS_start']], self.filtered_ecg_signal[indecies['QRS_start']], 'bx')
        plt.plot(self.t[indecies['P']], self.filtered_ecg_signal[indecies['P']], 'mo')
        plt.plot(self.t[indecies['T']], self.filtered_ecg_signal[indecies['T']], 'co')
        plt.plot(self.t[indecies['T_end']], self.filtered_ecg_signal[indecies['T_end']], 'mx')
        plt.plot(self.t, signal)
        #plt.figure(1)
        fig, ax1 = plt.subplots()
        ax1.plot(1e3 * (self.t[indecies['T_end']]) / self.fs, intervals['QT'], 'r')
        ax1.plot(1e3 * (self.t[indecies['T_end']]) / self.fs, intervals['QTP'], 'g')
        ax1.set_xlabel('time [s]')
        ax1.set_ylabel('QT')
        ax1.tick_params('y', colors='r')
        ax2 = ax1.twinx()
        ax2.plot(1e3 * (self.t[indecies['T_end']]) / self.fs, intervals['QTc'], 'b')
        ax2.plot(1e3 * (self.t[indecies['T_end']]) / self.fs, intervals['QTPc'], 'c')
        ax2.set_ylabel('QT_corrected')
        ax2.tick_params('y', colors='b')
        fig.tight_layout()
        plt.show(block=True)

    def _LT_QRS_detect(self, threshold_n, idx, peak_idx, peak_val, peak_indecies, zero_indecies, signal):
        threshold_n = 0.8 * threshold_n + 0.16 * abs(peak_val)
        peak_after = signal[peak_indecies[idx + 1]]
        peak_before = signal[peak_indecies[idx - 1]]

        zero_idx = bisect.bisect_right(zero_indecies, peak_idx)
        zero_idx -= 1 if abs(peak_before) > abs(peak_after) else 0

        R_idx = zero_indecies[zero_idx]
        S_idx = zero_indecies[zero_idx + 1]
        Q_idx = zero_indecies[zero_idx - 1]
        return Q_idx, R_idx, S_idx, threshold_n

    def _LT_RR_interval_calc(self, indecies, intervals, RR_interval_avg, R_idx):
        RR_interval = samples_to_ms(R_idx - indecies['R'][-1], self.fs)  # milliseconds
        intervals['RR'] = np.append(intervals['RR'], RR_interval)
        if RR_interval_avg is None:
            RR_interval_avg = RR_interval
        elif 1.5 * RR_interval_avg > RR_interval > 0.5 * RR_interval_avg:
            RR_interval_avg = 0.8 * RR_interval_avg + 0.2 * RR_interval
        return RR_interval_avg, RR_interval

    def _LT_QRS_start_detect(self, R_idx, Q_idx, indecies, peak_indecies, signal):
        if (R_idx - Q_idx) / self.fs < 0.08:  # 80 milliseconds
            indecies['Q'] = np.append(indecies['Q'], Q_idx)
            Q_i = peak_indecies[bisect.bisect_right(peak_indecies, Q_idx) - 1]
            threshold_q = signal[Q_i] / 2
            QRS_start_idx = Q_i
        else:
            R_i = peak_indecies[bisect.bisect_right(peak_indecies, R_idx) - 1]
            threshold_q = signal[R_i] / 5
            QRS_start_idx = R_i
        while abs(signal[QRS_start_idx]) > abs(threshold_q):
            QRS_start_idx -= 1
        return QRS_start_idx + 1

    def _LT_P_detect(self, indecies, QRS_start_idx, signal):
        bwind_idx = indecies['T_end'][-1]
        ewind_idx = QRS_start_idx
        window_max_idx = bwind_idx + np.argmax(signal[bwind_idx:ewind_idx])
        window_min_idx = bwind_idx + np.argmin(signal[bwind_idx:ewind_idx])
        P_idx = window_max_idx + np.argmax(self.filtered_ecg_signal[window_max_idx:window_min_idx])
        return P_idx

    def _LT_T_detect(self, RR_interval_avg, R_idx, indecies, zero_indecies, signal):
        bwind, ewind = (140, 500) if RR_interval_avg > 700 else (100, 0.7 * RR_interval_avg)
        bwind, ewind = ms_to_samples(bwind, self.fs), ms_to_samples(ewind, self.fs)
        bwind_idx = np.rint(R_idx + bwind).astype(int)  # window beginning
        ewind_idx = np.rint(R_idx + ewind).astype(int)  # window end
        window_max_idx = bwind_idx + np.argmax(signal[bwind_idx:ewind_idx])
        window_min_idx = bwind_idx + np.argmin(signal[bwind_idx:ewind_idx])

        if window_min_idx > window_max_idx:
            # If max is before min T is upward-downward or upward
            T_i = window_min_idx
            if abs(signal[window_max_idx]) > 4 * abs(signal[window_min_idx]):
                # If max is bigger than min -> T upward
                pass
            else:
                # If max is comparable to min -> T upward-downward
                pass
        else:
            # If min is before max T is downward, downward-upward or upward-downward
            # We're looking for min2 from max to the end of the window
            window_min2_idx = window_max_idx + np.argmin(signal[window_max_idx:ewind_idx])
            T_i = window_min2_idx
            if abs(signal[window_max_idx]) < 4 * abs(signal[window_min2_idx]):
                # If max is comparable to min2 -> T upward-downward
                pass
            else:
                # If max is bigger than min2
                if abs(signal[window_min_idx]) > 4 * abs(signal[window_max_idx]):
                    # If min is bigger than max -> T downward
                    pass
                else:
                    # If min is comparable to max -> T downward-upward
                    pass
        threshold_t = signal[T_i] / 2
        T_end_idx = T_i
        while T_end_idx < len(signal) and abs(signal[T_end_idx]) > abs(threshold_t):
            T_end_idx += 1
        T_idx = zero_indecies[bisect.bisect_right(zero_indecies, T_i) - 1]
        return (T_end_idx - 1), T_idx

    def _LT_intervals_calc(self, T_end_idx, T_idx, QRS_start_idx, RR_interval):
        QT_interval = samples_to_ms(T_end_idx - QRS_start_idx, self.fs)  # milliseconds
        QT_peak_interval = samples_to_ms(T_idx - QRS_start_idx, self.fs)  # milliseconds
        QT_corrected_interval = samples_to_ms((T_end_idx - QRS_start_idx) / np.sqrt(RR_interval), self.fs)
        QT_peak_corrected_interval = samples_to_ms((T_idx - QRS_start_idx) / np.sqrt(RR_interval), self.fs)
        return QT_interval, QT_peak_interval, QT_corrected_interval, QT_peak_corrected_interval

    @normalized
    def _lowpass(self, signal):
        b, a = butter(self.lowpass_order, self.lowpass_fc / self.nyquist_f)
        return filtfilt(b, a, signal)

    @normalized
    def _highpass(self, signal):
        b, a = butter(self.highpass_order, self.highpass_fc / self.nyquist_f, 'highpass')
        return filtfilt(b, a, signal)

    @normalized
    def _derivative(self, signal):
        return np.gradient(signal)

    @normalized
    def _squaring(self, signal):
        return signal**2

    @normalized
    def _moving_window_integral(self, signal):
        b = np.ones(self.integration_window_length + 1) / (self.integration_window_length + 1)
        a = [1]
        f = fft.rfftfreq(self.length, self.dt)
        w, gd = group_delay((b, a), f)
        shift = -np.rint(np.mean(gd)).astype(int)
        return np.roll(lfilter(b, a, signal), shift)

    def _spectrum(self, signal):
        return Spectrum(signal, self.t, self.length, self.dt)


def main():
    ecg = ECGAnalyzer()
    # for i in range(1, 10):
    #     ecg.highpass_order = i
    ecg.calculate()
    # plt.plot(1, 1)
    # plt.show(block=True)
    return
    m = (np.argmax(ecg.filtered_ecg_signal), max(ecg.filtered_ecg_signal))
    plt.figure(0)
    plt.plot(ecg.no_drift_ecg_signal)
    plt.show(block=False)
    plt.figure(1)
    plt.axhline(0, color='r')
    plt.plot(ecg.t, ecg.filtered_ecg_signal)
    plt.plot(ecg.t, ecg.no_drift_ecg_signal)
    plt.show(block=False)
    plt.figure(2)
    plt.plot(ecg.differentiated_ecg_signal)
    plt.show(block=False)
    plt.figure(3)
    plt.plot(ecg.squared_ecg_signal)
    plt.show(block=False)
    plt.figure(4)
    plt.plot(ecg.t, ecg.integrated_ecg_signal)
    plt.plot(ecg.t, ecg.filtered_ecg_signal)
    plt.show(block=True)


# ecg = ECGAnalyzer()
# ecg.calculate()
main()
