/*
Code to simulate the lock-in amplifier quantum measurements.
Authors: Andrii Sokolov, Elena Blokhina.
Equal1 Labs.
*/

#include <iostream>
#include <math.h>
#include <vector>
#include <random>
#include <numeric>
#include <fstream>
#include <omp.h>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

void pbar(float progress)
{
    int barWidth = 70;

    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i)
    {
        if (i < pos)
            std::cout << "=";
        else if (i == pos)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}

class LockInPulsing
{
private:
    const double total_time = 3.0;    // Total time [s]
    const double resolution = 20.0e6; // Resolution [1/s]
    const double f_lock_in = 10.0;    // Lock-in frequency [Hz]
    const double f_reset = 0.1e6;     // Reset frequency [Hz]
    const double f_rabi = 1.0e6;      // Rabi frequency [Hz]
    const double i_leak = 0.0e-12;    // Leakage current [A]
    const double i_rabi = 0.5e-12;    // Rabi frequency current [A]
    const double i_reset = 75.0e-12;  // Reset current [A]
    const double i_noise = 0.0;       // Noise current [A]
    const double tia_amp = 1e9;       // TIA amplification coefficient [V/A]
    const double Vref_amp = 1e-3;     // Amplitude of the reference voltage [V]
    const double mean_i_DUT = 0.0;    // Mean value of random current noise [A]
    const double std_i_DUT = 1e-11;   // STD of random current noise [A]
    double t_pulse;                   // Pulsing time [s]
    size_t number_of_points;          // Lengths of arrays
    std::vector<double> time;         // Array to hold time
    std::vector<double> Vref;         // Array to hold reference voltage
    std::vector<double> Vref2;        // Phase shifted reference voltage
    std::vector<double> iDUT;         // DUT current
    std::vector<double> VDUT;         // DUT voltage multiplied by reference voltage
    std::vector<double> VDUT2;        // DUT voltage multiplied by reference cos voltage
    double avg_X;
    double avg_Y;
    double avg_R;
    double avg_theta;

public:
    LockInPulsing(double t_p);
    ~LockInPulsing();
    void export_files(void);
    double Get_avg_X(void) { return (avg_X); };
    double Get_avg_Y(void) { return (avg_Y); };
    double Get_avg_R(void) { return (avg_R); };
    double Get_avg_theta(void) { return (avg_theta); };
};

int main(int argc, char **argv)
{
    size_t full_experiment = 1;
    if (full_experiment == 1)
    {
        omp_set_num_threads(20);
        size_t i, i_max, steps_completed;
        double t_pulse_max, t_iter;
        t_pulse_max = 5.0e-6;
        t_iter = 1.0e-7;
        i_max = floor(t_pulse_max / t_iter);
        std::ofstream file("export.csv");
        std::cout << "Calculating..." << std::endl
                  << "total number of iterations: \t" << i_max << std::endl;
#pragma omp parallel shared(file)
#pragma omp for
        for (i = 0; i < i_max; ++i)
        {
            LockInPulsing LiP(i * t_iter);
#pragma omp atomic
            ++steps_completed;

#pragma omp critical
            {
                file << i * t_iter << "," << LiP.Get_avg_X() << "," << LiP.Get_avg_Y() << "," << LiP.Get_avg_R() << "," << LiP.Get_avg_theta() << std::endl;
                pbar(steps_completed / i_max);
            }
        }
        std::cout << std::endl
                  << "Calculation finished" << std::endl;
    }
    else
    {
        LockInPulsing LiP(3e-6);
        LiP.export_files();
    }
    return (0);
}

LockInPulsing::LockInPulsing(double t_p)
{
    // Defining the parameters:
    t_pulse = t_p;
    number_of_points = floor(total_time * resolution);
    // Creating arrays:
    time.resize(number_of_points);
    Vref.resize(number_of_points);
    Vref2.resize(number_of_points);
    iDUT.resize(number_of_points);
    VDUT.resize(number_of_points);
    VDUT2.resize(number_of_points);
    // Creating the time array
    for (size_t i = 0; i < number_of_points; ++i)
    {
        time[i] = (float)i / resolution;
        Vref[i] = Vref_amp * cos(2.0 * M_PI * f_lock_in * time[i]);
        Vref2[i] = Vref_amp * cos(2.0 * M_PI * f_lock_in * time[i] + M_PI_2);
        if (sin(2.0 * M_PI * f_reset * time[i]) > 0)
        {
            iDUT[i] = i_reset;
        }
        else
        {
            iDUT[i] = i_leak;
        }
    }
    size_t i0, j;
    double t_pause = 1.0 / (2.0 * f_reset) - t_pulse;
    double t_total = 1.0 / (2.0 * f_reset);
    if (t_pulse > t_total)
    {
        t_pulse = t_total;
        t_pause = 0.0;
    }
    // Adding the rabi-oscillations
    for (size_t i = 0; i < number_of_points; ++i)
    {
        if ((sin(2.0 * M_PI * f_reset * time[i]) < 0) & (cos(2.0 * M_PI * f_lock_in * time[i]) > 0))
        {
            i0 = i;
            for (j = 0; time[j] < t_total; j += 1)
            {
                if (time[j] >= t_pause)
                    if (iDUT[i0 + j] != i_reset)
                        iDUT[i0 + j] += i_rabi * 0.5 * (cos(2.0 * M_PI * f_rabi * (time[j] - t_pause) - M_PI) + 1);
            }
            i = i0 + j;
        }
    }
    // Adding the random noise
    std::default_random_engine generator;
    std::normal_distribution<double> dist(mean_i_DUT, std_i_DUT);
    // Add Gaussian noise
    for (auto &x : iDUT)
    {
        x = x + dist(generator);
    }

    for (size_t i = 0; i < number_of_points; ++i)
    {
        VDUT[i] = iDUT[i] * tia_amp * Vref[i];
        VDUT2[i] = iDUT[i] * tia_amp * Vref2[i];
    }
    avg_X = std::accumulate(VDUT.begin(), VDUT.end(), 0.0) / VDUT.size();
    avg_Y = std::accumulate(VDUT2.begin(), VDUT2.end(), 0.0) / VDUT2.size();
    avg_R = sqrt(pow(avg_X, 2) + pow(avg_Y, 2));
    avg_theta = atan2(avg_Y, avg_X);
}

LockInPulsing::~LockInPulsing()
{
    time.clear();
    Vref.clear();
    Vref2.clear();
    iDUT.clear();
    VDUT.clear();
    VDUT2.clear();
}

void LockInPulsing::export_files(void)
{
    std::cout << "Exporting the csv file:" << std::endl;
    // Opening the export file:
    std::ofstream exportfile("export.csv");
    for (size_t i = 0; i < number_of_points; ++i)
    {
        // Exporting the file
        exportfile << time[i] << "," << Vref[i] << "," << Vref2[i] << "," << iDUT[i] << "," << VDUT[i] << std::endl;
        // Printing the progressbar
        if (i % 10000 == 0)
            pbar((double)i / (double)number_of_points);
    }
    std::cout << std::endl
              << "File successfully exported." << std::endl;
    exportfile.close();
}