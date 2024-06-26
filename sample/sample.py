from argparse import ArgumentParser

import emout
import matplotlib.pyplot as plt
from vdsolver.core import Particle, PhaseGrid

from vdsolverf.emses import get_probabilities


def parse_args():
    parser = ArgumentParser()

    parser.add_argument("--directory", "-d", default="./")
    parser.add_argument("--output", "-o", default="sample.png")

    return parser.parse_args()


def main():
    args = parse_args()

    data = emout.Emout(args.directory)

    NVX = 50
    NVZ = 50
    phase_grid = PhaseGrid(
        x=32,
        y=32,
        z=400,
        vx=(-50, 50, NVX),
        vy=0,
        vz=(-100, 100, NVZ),
    )

    phases = phase_grid.create_grid()
    particles = []
    for phase in phases.reshape(-1, phases.shape[-1]):
        pos = phase[:3].copy()
        vel = phase[3:].copy()
        pcl = Particle(pos, vel)
        particles.append(pcl)

    probs, ret_particles = get_probabilities(
        directory=args.directory,
        ispec=0,
        istep=-1,
        particles=particles,
        dt=data.inp.dt,
        max_step=100000,
        use_adaptive_dt=False,
    )

    phases = phases.reshape(NVZ, NVX, 6)
    VX = phases[:, :, 3]
    VZ = phases[:, :, 5]
    probs = probs.reshape(NVZ, -1)
    plt.pcolormesh(VX, VZ, probs, shading="auto")
    plt.colorbar()
    plt.gcf().savefig(args.output)


if __name__ == "__main__":
    main()
