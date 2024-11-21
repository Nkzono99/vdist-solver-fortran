import numpy as np


class Particle:
    """The Particle in the phase space."""

    def __init__(
        self, pos: np.ndarray, vel: np.ndarray, t: float = 0, periodic: bool = False
    ):
        """Initialize the Particle in the phase space.

        Parameters
        ----------
        pos : np.ndarray
            position
        vel : np.ndarray
            velocity
        t : float, optional
            time, by default 0
        periodic : bool, optional
            True when crossing the periodic boundary by default False

            This argment is used to visualize particle orbits (e.g. vdsolver.tools.plot.plot_periodic).
        """
        self.pos = pos
        self.vel = vel
        self.t = t
        self.periodic = periodic

    def __str__(self):
        return "Particle(p={pos}, v={vel}, t={t})".format(
            pos=self.pos,
            vel=self.vel,
            t=self.t,
        )

    @classmethod
    def create_prototype(cls, *args, **kwargs):
        pos = np.zeros(3)
        vel = np.zeros(3)
        return cls(pos, vel, *args, **kwargs)

    def craete_clone(self, pos: np.ndarray, vel: np.ndarray):
        return Particle(pos, vel)


class DustParticle:
    def __init__(
        self,
        charge: float,
        mass: float,
        radius: float,
        pos: np.ndarray,
        vel: np.ndarray,
    ):
        self.charge = charge
        self.mass = mass
        self.radius = radius
        self.pos = pos
        self.vel = vel
