# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import numpy as np
from PionScatLib import legP as legPOrig
# from matplotlib import pyplot as plt


def legP(x, n, deriv):
    if isinstance(x, np.ndarray):
        return np.array([legPOrig(a, n, deriv) for a in x])
    else:
        return legPOrig(x, n, deriv)


def main():
    for numDerivs in [1, 2]:
        testPass = True
        for n in range(9):

            def legPtmp(x):
                return legP(x, n=n, deriv=0)

            xs = np.arange(-1, 1, 0.1)
            ys1 = legP(xs, n=n, deriv=numDerivs)
            ys2 = deriv(legPtmp, order=numDerivs)(xs)

            # plt.scatter(xs, ys1)
            # plt.scatter(xs, ys2)

            for i, (y1, y2) in enumerate(zip(ys1, ys2)):
                diff = abs(y1 - y2)
                if diff > 1e-4:
                    round = np.round(xs[i], 4)
                    print(
                        f"Issue at i={i:2d}, x={round:+10.6f}, y1={y1:+10.6f}, y2={y2:+10.6f}, n={n:2d},   diff={diff}"
                    )
                    testPass = False

        # plt.show()
        print(f"For numDerivs={numDerivs}, Test Pass=", testPass)


def deriv(foo, order=1):
    h = 0.00001

    def fooPrime(x):
        return (foo(x + h) - foo(x - h)) / (2 * h)

    order = order - 1
    if order == 0:
        return fooPrime
    else:
        return deriv(fooPrime, order)


if __name__ == "__main__":
    main()
