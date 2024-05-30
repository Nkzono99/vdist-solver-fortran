from argparse import ArgumentParser

import emout
import matplotlib.pyplot as plt
from vdsolver.core import Particle

from vdsolverf import get_backtrace


def parse_args():
    parser = ArgumentParser()

    parser.add_argument("--directory", "-d", default="./")
    parser.add_argument("--output", "-o", default="backtrace-sample.png")

    return parser.parse_args()


def main():
    args = parse_args()

    data = emout.Emout(args.directory)

    position = [32, 32, 400]
    velocity = [0, 0, -10]
    max_step = 100000

    particle = Particle(position, velocity)

    ts, positions, velocities = get_backtrace(
        directory=args.directory,
        ispec=0,
        istep=-1,
        particle=particle,
        dt=data.inp.dt,
        max_step=max_step,
        use_adaptive_dt=False,
    )

    data.phisp[-1, :, 32, :].plot(use_si=False)

    plt.plot(positions[:, 0], positions[:, 2])
    plt.gcf().savefig(args.output)


if __name__ == "__main__":
    main()
