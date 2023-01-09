"""
The backend module that provides the simulator of the lock-in measurements
Authors: Andrii Sokolov, Elena Blokhina
"""

import matplotlib.pyplot as plt
import numpy as np
from console_progressbar import ProgressBar
from scipy import signal


class Lockin:
    """
    This class simulates the Lock-in measurements of a DUT with given i-V curve and
    non-stationary noise.
    """

    def __init__(
        self,
        n_time: int = 20000,
        n_rc: int = 30,
        f_ref: float = 3.14,
        a_ref: float = 0.001,
        v_bias: float = 0.0,
        pref: float = 0.0,
        rc_time: float = 1.5,
        n_filter: float = 2,
        noise_1: float = 0.2,
        noise_2: float = 0.5,
        noise_3: float = 0.5,
        non_stat_noise: float = 0.0,
        n_telegraph: float = 100,
        sr570_gain: float = 1e-9,
        sr570_f_file: str = "SR570_Gain_f.csv",
        sr865a_sens: float = 100,
    ):
        # Number of points in time-series:
        self.n_time = n_time
        # Number of time-constant periods:
        self.n_rc = n_rc
        # Frequency of the reference signal [Hz]:
        self.f_ref = f_ref
        # Amplitude of the reference signal [V]:
        self.amp_ref = a_ref
        # Phase of the reference signal [deg]:
        self.p_ref = pref * np.pi / 180
        # Time constant [s]:
        self.rc_time = rc_time
        # Cut-off frequency [Hz]:
        self.f_cut_off = 1.0 / (2.0 * np.pi * self.rc_time)
        # Filter order:
        self.n_filter = n_filter
        # Array of time:
        self.time = np.linspace(0, self.n_rc * self.rc_time, self.n_time)
        # Reference signal (Y):
        self.sin_ref = v_bias + self.amp_ref * np.sin(
            2.0 * np.pi * self.f_ref * self.time + self.p_ref
        )
        # Reference signal (X):
        self.cos_ref = v_bias + self.amp_ref * np.cos(
            2.0 * np.pi * self.f_ref * self.time + self.p_ref
        )
        # Adding the noise:
        self.sin_ref += np.array(
            [
                self.amp_ref * noise_1 * (2 * np.random.random() - 1)
                for a in range(self.n_time)
            ]
        )
        self.cos_ref += np.array(
            [
                self.amp_ref * noise_1 * (2 * np.random.random() - 1)
                for a in range(self.n_time)
            ]
        )
        # Relative noise on the pre-amp:
        self.noise_2 = noise_2
        # Relative noise on the lock-in amp:
        self.noise_3 = noise_3
        # Non-stationary noise coefficient:
        self.non_stat_noise = non_stat_noise
        # Non-stationary noise number of transitions:
        self.n_telegraph = n_telegraph
        # Telegraph_noise array:
        self.telegraph = np.zeros(self.n_time)
        # TiA SR570A Gain:
        self.sr570_gain = sr570_gain
        # Reading the frequency characteristic of the TiA:
        [self.f_sr570, self.g_sr570] = np.genfromtxt(sr570_f_file, delimiter=",").T
        # The current after DUT:
        self.i_dut = np.array([])
        # The voltage after TiA:
        self.u_tia = np.array([])
        # The voltage after lock-in amp:
        self.u_lockin = np.array([])
        # The sensitivity of the SR865a:
        self.sr865a_sens = sr865a_sens
        # The outputs of the lockin:
        self.sr865a_x = np.array([])
        self.sr865a_y = np.array([])
        self.sr865a_r = np.array([])
        self.sr865a_theta = np.array([])

    def sr570_tia(self, current: float):
        """Function that simulates the behavior of
        the SRS SR570 Trans-impedance Amplifier

        Args:
            current (float): the input current

        Returns:
            (float): output voltage
        """
        return current / (
            self.sr570_gain
            / 10
            ** np.where(
                self.f_ref < self.f_sr570[0],
                self.g_sr570[0],
                np.where(
                    self.f_ref > self.f_sr570[-1],
                    self.g_sr570[-1],
                    np.interp(self.f_ref, self.f_sr570, self.g_sr570),
                ),
            )
        )

    def lock_in_measurement(self, dut_function: object):
        """Function that simulates the lockin measurement
        adding the amplification of TiA and lockin amplifier

        Args:
            dut_function (object): Function that simulates DUT
        """
        # Calculating the reference signal after the DUT:
        self.i_dut = dut_function(self.sin_ref)
        # Define maximum and minimum of the current:
        i_min = np.min(self.i_dut)
        i_max = np.max(self.i_dut)
        # Creating the random telegraph time-instances:
        tsw = np.sort(np.random.rand(self.n_telegraph)) * self.time[-1]
        # Switch:
        sw = True
        # index of telegraph
        j = 0
        for k in range(self.n_time):
            if self.time[k] > tsw[j]:
                if j < len(tsw) - 1:
                    j += 1
                sw = not sw
            if sw:
                self.telegraph[k] = self.non_stat_noise * i_max
            else:
                self.telegraph[k] = self.non_stat_noise * i_min
        self.i_dut += self.telegraph
        # Calculating the voltage after the TiA:
        self.u_tia = self.sr570_tia(self.i_dut)
        # Adding the noise to the TiA signal:
        u_tia_max = np.max(self.u_tia)
        u_tia_min = np.min(self.u_tia)
        self.u_tia += np.array(
            [
                self.noise_2 * (u_tia_max - u_tia_min) * (2 * np.random.random() - 1)
                for a in range(self.n_time)
            ]
        )
        # Applying the amplifier of the Lockin:
        self.u_lockin = (self.u_tia / self.sr865a_sens) * 10
        # Applying the noise to the Lockin:
        u_lockin_max = np.max(self.u_lockin)
        u_lockin_min = np.min(self.u_lockin)
        self.u_lockin += np.array(
            [
                self.noise_3
                * (u_lockin_max - u_lockin_min)
                * (2 * np.random.random() - 1)
                for a in range(self.n_time)
            ]
        )
        # Removing the DC-components:
        no_dc_sig = self.u_lockin - np.mean(self.u_lockin)
        no_dc_sin_ref = self.sin_ref - np.mean(self.sin_ref)
        no_dc_cos_ref = self.cos_ref - np.mean(self.cos_ref)
        # Mixing the signal:
        mix_x = np.multiply(no_dc_sig, no_dc_sin_ref)
        mix_y = np.multiply(no_dc_sig, no_dc_cos_ref)
        # Setting up the filter:
        sos = signal.butter(
            self.n_filter,
            self.f_cut_off,
            "low",
            output="sos",
            fs=1 / (self.time[1] - self.time[0]),
        )
        # Applying the lowpass filter:
        # X-output
        self.sr865a_x = signal.sosfilt(sos, mix_x)
        # Y-output
        self.sr865a_y = signal.sosfilt(sos, mix_y)
        # R-output
        self.sr865a_r = np.sqrt(np.square(self.sr865a_x) + np.square(self.sr865a_y))
        # Theta-output
        self.sr865a_theta = np.arctan2(np.real(self.sr865a_x), np.real(self.sr865a_y))


def dut_resistor(v_ds: float, resistance: float = 1e12):
    """DUT function of the simple resistor

    Args:
        v_ds (float): input voltage on the device
        r (float, optional): resistance. Defaults to 1e12.

    Returns:
        float: output current
    """
    return v_ds / resistance


class DUTSemiconductor:
    """Class that simulates the behaviour of the semiconductor DUT"""

    def __init__(
        self,
        vg_max: list = (0.0, 5e-3),
        g_max: list = (2e-9, 1e-9),
        w_max: list = (1e-6, 2e-6),
        exp_a: float = 6e-11,
        exp_b: float = 30.0,
        n_points: int = 10001,
        vds_min: float = -2.5e-2,
        vds_max: float = 2.5e-2,
    ):
        # Initial current
        self.current = np.zeros(n_points)
        # Generating the conductance versus v_DS
        self.v_ds = np.linspace(vds_min, vds_max, n_points)
        self.g_ds = exp_a * exp_b * np.exp(exp_b * self.v_ds)
        for vg, g, w in zip(vg_max, g_max, w_max):
            self.g_ds += g * np.exp(-np.square(self.v_ds - vg) / w)
        # integrating the current array:
        for iter_vds in range(len(self.v_ds)):
            self.current[iter_vds] = np.sum(self.g_ds[:iter_vds]) * (
                self.v_ds[1] - self.v_ds[0]
            )

    def g_interp(self, v_in: float):
        """interpolation function that returnes the conductancre

        Args:
            v_in (float): input voltage

        Returns:
            float: output conductance
        """
        return np.interp(v_in, self.v_ds, self.g_ds)

    def i_interp(self, v_in: float):
        """interpolation function that returnes the conductancre

        Args:
            v_in (float): input voltage

        Returns:
            float: output current
        """
        return np.interp(v_in, self.v_ds, self.current)


def telegraphic_noise_coefficient(v_in: float, vb_max: list = (5e-3,), w_max=(3e-5,)):
    """Telegraphic noise coefficients returns coefficient from 0 to 1 depending on the bias voltage

    Args:
        v_in (float): input voltage
        vb_max (list, optional): maximums. Defaults to (-1e-2, -5e-3).
        w_max (tuple, optional): width of the maximums. Defaults to (3e-6, 3e-6).

    Returns:
        float: from 0 to 1
    """
    out = 0.0
    for v1, wb in zip(vb_max, w_max):
        out += np.exp(-np.square(v_in - v1) / wb)
    return out


if __name__ == "__main__":
    pb = ProgressBar(
        total=100,
        prefix="Progress",
        decimals=3,
        length=50,
        fill="#",
        zfill="-",
    )
    SC = DUTSemiconductor()
    r_out = []
    t_out = []
    n_vbiases = 100
    v_bias_test = np.linspace(-0.02, 0.02, n_vbiases)
    for vb, i in zip(v_bias_test, range(n_vbiases)):
        Li = Lockin(v_bias=vb)
        Li.lock_in_measurement(SC.i_interp)
        r_out.append(Li.sr865a_r[-1])
        t_out.append(Li.sr865a_theta[-1])
        pb.print_progress_bar(i)
    plt.plot(v_bias_test, r_out, "k.")
    ax = plt.twinx()
    ax.plot(v_bias_test, t_out, "r", linewidth=0.5)
    plt.show()
