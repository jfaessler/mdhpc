import numpy as np
import matplotlib.pyplot as plt


# def f(x, a, b):
#     return x ** a + b


x = [
    10,
    50,
    100,
    150,
    200,
    500
]
y = [
    0.157,
    2.576,
    7.040,
    22.526,
    40.609,
    60 * 3 + 16.841
]

plt.scatter(x, y)
plt.show()
