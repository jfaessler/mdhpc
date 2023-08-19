import glob

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.fft import fft, fftfreq

from typing import List


class Run:
    def __init__(self, file_name: str):
        self.file_name = file_name
        with open(file_name) as file:
            # self.comment_line = file.readline().replace('\n', '')

            self.params = {}
            self.data = {}
            label_line = None

            for line in file.readlines():
                if line[0] == '#':
                    if line[0:8] == "#PARAMS:":
                        params_line = line[8:].replace('\n', '').split(',')
                        for param in params_line:
                            k, v = param.split('=')
                            self.params[k] = float(v)
                    continue
                if label_line is None:
                    label_line = line.replace('\n', '').split(',')
                    for label in label_line:
                        self.data[label] = []
                    continue

                split = line.replace('\n', '').split(',')
                if len(split) != len(self.data.keys()):
                    print("Wrong number of delimiters in line: " + str(line))
                    line = file.readline()
                    continue
                for i, value in enumerate(self.data.values()):
                    value.append(float(split[i]))


def energy_time(report: Run):
    plt.plot(report.data['Step'], report.data['Temperature'])
    plt.show()


def temperature_step(report: Run, step_scale=50):
    x = np.array(report.data['Step'])
    # plt.plot(x, report.data['Average Temperature'])
    # plt.show()
    kinetic = report.data['Kinetic']
    plt.plot(x, kinetic)
    # plt.show()
    pot = report.data['Potential']
    total = []
    for i in range(len(pot)):
        pot[i] += 1050  # Potential is shifted for visualization
        total.append(pot[i] + kinetic[i])
    plt.plot(x, pot)
    plt.plot(x, total)
    plt.show()

    def f(x, a, b):
        return a * x + b

    popt, pcov = curve_fit(f, x, total)
    plt.plot(x, f(x, *popt), "b+")
    plt.plot(x, total, "r-")
    plt.show()
    print(popt)

    # plt.plot(report.data['Average Temperature'], total)
    # plt.plot(total, report.data['Average Temperature'])
    # plt.show()


def cap_size(reports):
    # TODO check that system is fully equalized
    def linear(x, a, b):
        return a * x + b

    def heat_curve(x, t1, t2, s1, s2, height):
        return np.piecewise(x, [x < t1, (t1 <= x) & (x <= t2), t2 < x],
                            [lambda x: s1 * x + height,
                             lambda x: s1 * t1 + height + 0 * x,
                             lambda x: s1 * t1 + height + s2 * (x - t2)])

    size = []
    heat_capacity = []
    melting_point = []
    latent_heat = []
    for report in reports:
        q = np.array(report.data['Step']) * report.params['delta_q']
        T = np.array(report.data['Average Temperature'])
        size.append(report.params['size'])
        heat_opt, pcov = curve_fit(heat_curve, q, T, p0=[75000 + size[-1] * 7, 100000 + size[-1] * 7, 1, 1, 1])
        latent_heat.append(heat_opt[1] - heat_opt[0])
        plt.plot(q, heat_curve(q, *heat_opt), label='Fitted Trajectory')
        plt.plot(q, T, label='Observed Temperature')
        plt.legend()
        plt.xlabel('Energy Added (eV)')
        plt.ylabel('Temperature (K)')
        heat_capacity.append(heat_opt[2])
        melting_point.append(heat_curve(heat_opt[0], *heat_opt))
        plt.show()
    plt.plot(size, heat_capacity)
    plt.xlabel("Cluster Size N (# atoms)")
    plt.ylabel("Heat Capacity C (eV / K)")
    plt.show()
    plt.plot(size, latent_heat)
    plt.xlabel("Cluster Size N (# atoms)")
    plt.ylabel("Latent Heat (eV)")
    plt.show()
    plt.plot(size, melting_point)
    plt.xlabel("Cluster Size N (# atoms)")
    plt.ylabel("Melting Point (K)")
    plt.show()


def stress_strain(report):
    s = np.array(report.data['Step'])
    s /= 1000
    plt.plot(s, report.data['Stress'])
    plt.show()
    plt.plot(report.data['Strain'], report.data['Stress'])
    plt.show()


def energy_time_4(reports: List[Run], max_timestep: float):
    for report in reports:
        if report.params['Timestep'] < max_timestep:
            plt.plot(report.data['Time'], report.data['Total Energy'],
                     label=('Timestep: ' + str(report.params['Timestep'])))
    plt.legend()
    plt.xlabel('Time')
    plt.ylabel('Total Energy')
    plt.tight_layout()
    plt.show()


def scaling_neighbor():
    # TODO outlier
    def quadratic(x, a):
        return a * (x**2)
    # Hardcode values from scaling.sh run
    cluster_size = np.array([4 ** 3, 5 ** 3, 6 ** 3, 7 ** 3, 8 ** 3, 9 ** 3])
    time5 = np.array([5.508, 20.815, 60 + 1.818, 2 * 60 + 36.214, 5 * 60 + 46.151, 13 * 60 + 17.415])
    # cluster_size = np.array([4 ** 3, 5 ** 3, 6 ** 3, 7 ** 3, 8 ** 3])
    # time5 = np.array([5.508, 20.815, 60 + 1.818, 2 * 60 + 36.214, 5 * 60 + 46.151])
    time6 = np.array([6.664, 24.605, 60 + 3.911, 60 * 2 + 19.504, 3 * 57.639, 6 * 60 + 56.023])
    plt.plot(cluster_size, time5)
    plt.plot(cluster_size, time6)
    plt.xlabel("Cluster Size (# atoms)")
    plt.ylabel("Runtime (s)")
    plt.show()
    popt, pcov = curve_fit(quadratic, cluster_size, time5)
    plt.loglog(cluster_size, time5)
    # plt.loglog(cluster_size, quadratic(cluster_size, *popt))
    plt.loglog(cluster_size, time6)
    plt.xlabel("Cluster Size (# atoms)")
    plt.ylabel("Runtime (s)")
    plt.show()


if __name__ == '__main__':
    scaling_neighbor()
    files = []
    reports = []
    for filename in glob.glob("*.csv"):
        files.append(filename)
    files.sort()
    for filename in files:
        reports.append(Run(filename))
    # for report in reports:
    # energy_time(report)
    # temperature_step(report)
    # stress_strain(report)
    cap_size(reports)
    files = []
    reports = []
    for filename in glob.glob("*.timestep"):
        files.append(filename)
    files.sort()
    for filename in files:
        reports.append(Run(filename))
    # energy_time_4(reports, 0.03)
    # energy_time_4(reports, 0.05)
    # energy_time_4(reports, 1.0)
