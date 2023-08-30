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
    plt.tight_layout()
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
    plt.tight_layout()
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
        # Rescale energy to be 0 at the start
        initial_energy = report.data['Kinetic'][0] + report.data['Potential'][0]
        e_added = np.array(np.array(report.data['Kinetic']) + np.array(report.data['Potential']) - initial_energy)
        T = np.array(report.data['Average Temperature'])
        plt.plot(T, report.data['Potential'], label='Potential Energy')
        plt.plot(T, np.array(report.data['Kinetic']) + np.array(report.data['Potential']), label='Total Energy')
        plt.xlabel('Temperature (K)')
        plt.ylabel('Energy (eV)')
        plt.legend()
        plt.tight_layout()
        plt.show()
        size.append(report.params['size'])
        heat_opt, pcov = curve_fit(heat_curve, e_added, T, p0=[size[-1] / 10 - 20, size[-1] / 5 - 20, 1, 1, 1])
        latent_heat.append(heat_opt[1] - heat_opt[0])
        plt.plot(e_added, heat_curve(e_added, *heat_opt), label='Fitted Equation')
        plt.plot(e_added, T, label='Observed Temperature')
        plt.legend()
        plt.xlabel('Energy Added (eV)')
        plt.ylabel('Temperature (K)')
        heat_capacity.append(1 / heat_opt[2])  # Approximate dE/dT with the fitted slope (which is dT/dE)
        melting_point.append(heat_curve(heat_opt[0], *heat_opt))
        plt.tight_layout()
        plt.show()

    size = np.array(size)

    plt.plot(size, heat_capacity, label="Heat Capacity")
    plt.xlabel("Cluster Size N (# atoms)")
    plt.ylabel("Heat Capacity C (eV / K)")
    popt, pcov = curve_fit(linear, size, heat_capacity)
    plt.plot(size, linear(size, *popt), 'k--', label='Fitted')
    plt.legend()
    plt.show()
    mass = size * 20413.15887 / 103.6
    popt, pcov = curve_fit(linear, mass, heat_capacity)
    plt.plot(mass, linear(mass, *popt), 'k--', label='Fitted')
    print("Specific Heat Capacity: {} eV/g/molK".format(popt[0]))
    plt.tight_layout()
    plt.show()

    plt.plot(size, latent_heat)
    plt.xlabel("Cluster Size N (# atoms)")
    plt.ylabel("Latent Heat (eV)")
    plt.tight_layout()
    plt.show()

    plt.plot(size, melting_point)
    plt.xlabel("Cluster Size N (# atoms)")
    plt.ylabel("Melting Point (K)")
    plt.tight_layout()
    plt.show()


def stress_strain(report):
    s = np.array(report.data['Step'])
    s /= 1000
    plt.plot(s, report.data['Stress'])
    plt.xlabel("Simulation Step")
    plt.ylabel("Stress (eV / \u212b^3)")
    plt.tight_layout()
    plt.show()
    plt.plot(report.data['Strain'][6:], report.data['Stress'][6:])
    plt.xlabel("Strain \u03b5")
    plt.ylabel("Stress \u03c3 (eV / \u212b^3)")
    plt.tight_layout()
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
    def quadratic(x, a):
        return a * (x ** 2)

    def linear(x, a):
        return a * x

    # Hardcode values from scaling.sh run
    cluster_size = np.array([4 ** 3, 5 ** 3, 6 ** 3, 7 ** 3, 8 ** 3, 9 ** 3, 10 ** 3, 11 ** 3, 12 ** 3])
    time5 = np.array([5.798, 22.119, 60 + 6.095, 2 * 60 + 47.051, 5 * 60 + 51.451, 12 * 60 + 7.075, 22 * 60 + 44.844,
                      40 * 60 + 1.476, 67 * 60 + 46.742])
    time6 = np.array([6.775, 26.213, 60 + 6.307, 2 * 60 + 15.037, 3 * 60 + 55.054, 6 * 60 + 55.267, 10 * 60 + 10.301,
                      15 * 60 + 46.771, 21 * 60 + 26.981])

    plt.plot(cluster_size, time5, label='Without Neighbor Lists')
    plt.plot(cluster_size, time6, label='With Neighbor Lists')
    plt.xlabel("Cluster Size (# atoms)")
    plt.ylabel("Runtime (s)")
    plt.legend()
    plt.tight_layout()
    plt.show()

    popt, pcov = curve_fit(quadratic, cluster_size, time5)
    plt.loglog(cluster_size, quadratic(cluster_size, *popt), 'k--', label='O(N^2) Scaling')
    plt.loglog(cluster_size, time5, label='Without Neighbor Lists')

    # popt, pcov = curve_fit(linear, cluster_size, time6)
    # plt.loglog(cluster_size, linear(cluster_size, *popt), 'b--', label='O(N) Scaling')
    plt.loglog(cluster_size, time6, label='With Neighbor List')

    plt.xlabel("Cluster Size (# atoms)")
    plt.ylabel("Runtime (s)")
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    # scaling_neighbor()
    files = []
    reports = []
    for filename in glob.glob("*.csv"):
        files.append(filename)
    files.sort()
    for filename in files:
        reports.append(Run(filename))
    for report in reports:
        # energy_time(report)
        # temperature_step(report)
        stress_strain(report)
    # cap_size(reports)
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
