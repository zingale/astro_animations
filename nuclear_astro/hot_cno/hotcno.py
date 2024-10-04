#/usr/bin/env python

import matplotlib
matplotlib.use('agg')

import matplotlib.pyplot as plt
import numpy as np

import pynucastro as pyna

def main():

    rl = pyna.ReacLibLibrary()

    rate_names = ["c12(p,g)n13",
                  "c13(p,g)n14",
                  "n13(,)c13",
                  "n13(p,g)o14",
                  "n14(p,g)o15",
                  "n15(p,a)c12",
                  "o14(,)n14",
                  "o15(,)n15",
                  "o14(a,p)f17",
                  "f17(p,g)ne18",
                  "ne18(,)f18",
                  "f18(p,a)o15"]

    rates = rl.get_rate_by_name(rate_names)

    rc = pyna.RateCollection(rates=rates)

    comp = pyna.Composition(rc.get_nuclei())
    comp.set_solar_like()

    rho = 100

    Ts = np.logspace(7, np.log10(8.e8), 500)

    dpi = 150

    for n, T in enumerate(Ts):
        fig = rc.plot(rho=rho, T=T, comp=comp, dpi=dpi,
                      ydot_cutoff_value=1.e-30,
                      show_small_ydot=True,
                      size=(1280, 720))

        ax = fig.axes[0]

        T_exp = int(np.log10(T))
        T_sig = T / 10**T_exp

        ax.text(0.05, 0.9, rf"$T = {T_sig:5.2f} \times 10^{T_exp}$ K",
                transform=ax.transAxes)
        ax.text(0.95, 0.04, "https://github.com/pynucastro/",
                horizontalalignment="right",
                transform=fig.transFigure)

        fig.savefig(f"hotcno_{n:03d}.png", dpi=dpi)
        plt.close(fig)

if __name__ == "__main__":
    main()
