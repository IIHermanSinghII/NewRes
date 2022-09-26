import sys
import unittest
from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from asap3 import Trajectory
from ase import units

from ase.calculators.emt import EMT
from md import calcenergy


class MdTests(unittest.TestCase):
    def test_calcenergy(self):
        size = 3
        atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                                  symbol='Cu',
                                  size=(size, size, size),
                                  pbc=True)
        atoms.calc = EMT()
        MaxwellBoltzmannDistribution(atoms, temperature_K=300)

        # We want to run MD with constant energy using the VelocityVerlet algorithm.
        dyn = VelocityVerlet(atoms, 5 * units.fs)  # 5 fs time step.
        traj = Trajectory("cu_traj", "w", atoms)
        dyn.attach(traj.write, interval=10)

        for i in range(10):
            dyn.run(3)
            epot, ekin, ins_temp, e_tot = calcenergy(atoms)
            print(epot)

        self.assertGreater(epot, 0)


if(__name__ == '__main__'):
    tests = [unittest.TestLoader().loadTestsFromTestCase(MdTests)]
    testsuite = unittest.TestSuite(tests)
    result = unittest.TextTestRunner(verbosity=0).run(testsuite)
    sys.exit(not result.wasSuccessful())
