"""
########################################################
#
# Model of the Lock-in measurement test for the pulsing
# experiment.
# A Sokolov E Blokhina
# Equal1 labs
#
########################################################
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import butter, sosfilt, square


class LockInPulsing:
    """
    This class simulates the Lock-in measurements of the pulsing experiment
    """

    def __init__(
        self,
        t_tot: float = 3.0,
        resolution: int = 10e6,
        f_lock_in: float = 10.0,
        v_amp_lock_in: float = 1e-3,
        f_reset: float = 0.1e6,
        f_rabi: float = 1.0e6,
        sigma_min: float = 0,
        sigma_amp: float = 40e-12,
        sigma_shift: float = 1e-6,
        t_pulse: float = 3e-6,
        tia_coefficient: float = 1e-3 / 1e-15,
    ):
        # Total time [s]
        self.t_tot = t_tot
        # Time resolution [points/s]
        self.resolution = resolution
        # Total number of points
        self.n_tot = int(t_tot * resolution)
        # Time array
        self.t_array = np.linspace(0, self.t_tot, self.n_tot)
        # Lock-in frequency [Hz]
        self.f_lock_in = f_lock_in
        # Lock-in reference amplitude [V]
        self.v_amp_lock_in = v_amp_lock_in
        # Lock-in reference signal [V]
        self.v_ref_cos = np.where(
            np.cos(2 * np.pi * self.f_lock_in * self.t_array) >= 0,
            self.v_amp_lock_in,
            -self.v_amp_lock_in,
        )
        self.v_ref_sin = np.where(
            np.sin(2 * np.pi * self.f_lock_in * self.t_array) >= 0,
            self.v_amp_lock_in,
            -self.v_amp_lock_in,
        )
        # Reset frequency [Hz]
        self.f_reset = f_reset
        # Rabi frequency [Hz]
        self.f_rabi = f_rabi

        # Minimal conductance [Si]
        self.sigma_min = sigma_min
        # Conductance variation [Si]
        self.sigma_amp = sigma_amp
        # Pulse time [s]
        self.t_pulse = t_pulse
        # Pulse time shift [s]
        self.sigma_shift = sigma_shift

        # Conductance model [Si]
        self.sigma = np.where(
            np.cos(2 * np.pi * self.f_lock_in * self.t_array) >= 0.0,
            self.sigma_min
            + self.sigma_amp
            * (
                np.where(
                    square(
                        2.0 * np.pi * self.f_reset * (self.t_array - sigma_shift),
                        duty=t_pulse * f_reset,
                    )
                    > 0.0,
                    np.cos(
                        2.0 * np.pi * self.f_rabi * (self.t_array - sigma_shift) - np.pi
                    ),
                    -1.0,
                )
                + 1
            )
            / 2,
            self.sigma_min,
        )
        self.sigma += self.sigma_amp * (
            (
                square(
                    2.0 * np.pi * self.f_reset * (self.t_array + sigma_shift) + np.pi,
                    duty=0.5,
                )
                + 1
            )
            / 2
        )
        # Current after DUT [A]
        self.iDUT = self.v_ref_cos * self.sigma
        # TiA coefficient [V/A]
        self.tia_coeff = tia_coefficient
        # Voliage in Lock-in [V]
        self.v_in = self.iDUT * self.tia_coeff
        # Multiplications [V]
        self.cos_multi = np.multiply(self.v_ref_cos, self.v_in)
        self.sin_multi = np.multiply(self.v_ref_sin, self.v_in)
        # Export X, Y, R, Theta
        self.x_out = 0.0
        self.y_out = 0.0
        self.r_out = 0.0
        self.theta_out = 0.0

    def plot_signals(self):
        """
        Method displays all initialized signals
        """
        _, axs = plt.subplots(nrows=5, figsize=(16, 10))
        axs[0].plot(self.t_array, 1e3 * self.v_ref_cos)
        axs[0].plot(self.t_array, 1e3 * self.v_ref_sin)
        axs[0].set_ylabel("$V_{ref}$, mV")
        axs[1].plot(self.t_array, self.sigma * 1e9, ".--")
        axs[1].set_xlim([0, 5 / self.f_reset])
        axs[1].set_ylabel("$\\sigma$, nSi")
        axs[2].plot(self.t_array, self.iDUT * 1e12, ".--")
        axs[2].set_xlim([0, 5 / self.f_reset])
        axs[2].set_ylabel("$i_{DUT}$, fA")
        axs[3].plot(self.t_array, self.v_in)
        axs[3].set_xlim([0, 5 / self.f_reset])
        axs[4].plot(self.t_array, self.cos_multi)
        axs[4].plot(self.t_array, self.sin_multi)
        for i in range(5):
            axs[i].set_xlabel("t, s")
        plt.tight_layout()

        plt.show()

    def ideal_filter(self):
        """
        Method returns average of signals
        X
        Y
        R
        Theta
        """
        self.x_out = np.mean(self.cos_multi)
        self.y_out = np.mean(self.sin_multi)
        self.r_out = np.sqrt(self.x_out**2 + self.y_out**2)
        self.theta_out = np.arctan2(self.x_out, self.y_out)
        return [self.x_out, self.y_out, self.r_out, self.theta_out]


if __name__ == "__main__":
    # LP = LockInPulsing(t_pulse=3e-6)
    # LP.plot_signals()
    tpulse = np.linspace(0, 3e-6, 301)
    print(tpulse)
    x = np.zeros(301)
    y = np.zeros(301)
    r = np.zeros(301)
    th = np.zeros(301)
    for i, tp in enumerate(tpulse):
        LP = LockInPulsing(t_pulse=tp)
        print(i, tp, LP.ideal_filter())
        x[i] = LP.x_out
        y[i] = LP.y_out
        r[i] = LP.r_out
        th[i] = LP.theta_out

    plt.plot(tpulse, x)
    plt.plot(tpulse, y)
    plt.plot(tpulse, r)
    plt.savefig("plot_amplitudes_3.png")
    plt.show()
    plt.plot(tpulse, th)
    plt.savefig("plot_phase_3.png")
    plt.show()
