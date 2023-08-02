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
    def linear(x, a, b):
        return a * x + b

    def heat_curve(x, t1, t2, s1, s2, height):
        return np.piecewise(x, [x < t1, (t1 <= x) & (x <= t2), t2 < x],
                            [lambda x: s1 * x + height,
                             lambda x: s1 * t1 + height + 0 * x,
                             lambda x: s1 * t1 + height + s2 * (x - t2)])

    def cap_fn(x, cap, height):
        return cap * x + height

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
        plt.plot(q, heat_curve(q, *heat_opt))
        plt.plot(q, T)
        plt.show()
        kinetic = np.array(report.data['Kinetic'])
        pot = np.array(report.data['Potential'])
        total = kinetic + pot
        cap_opt, pcov = curve_fit(cap_fn, q, total)
        # plt.plot(q, total, "r+")
        # plt.plot(q, cap_fn(q, *cap_opt))
        # plt.show()
        heat_capacity.append(cap_opt[0])
        melting_point.append(heat_curve(heat_opt[0], *heat_opt))
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
    plt.plot(report.data['Step'], report.data['Stress'])
    plt.show()
    plt.plot(report.data['Strain'], report.data['Stress'])
    plt.show()


if __name__ == '__main__':
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
