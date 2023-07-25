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
            params_line = None
            label_line = None

            for line in file.readlines():
                if line[0] == '#':
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

def temperature_step(report: Run, step_scale = 50):
    y = report.data['Temperature']
    x = np.arange(len(y)) * step_scale
    plt.plot(x, y)
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
        temperature_step(report)
