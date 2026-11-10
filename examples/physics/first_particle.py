"""Calculate the energy and force for one particle in a harmonic potential.

This is the smallest complete Moleculab example.  It uses the parameters from
the first milestone tutorial and prints the same quantities that are derived by
hand there:

    K = 0.06 kJ/mol
    U = 0.50 kJ/mol
    E = 0.56 kJ/mol
    F = (-10, 0, 0) kJ/(mol nm)
    a = (-0.8333, 0, 0) nm/ps^2

Run it with:

    uv run python examples/physics/first_particle.py
"""

from moleculab import (
    HarmonicPotential,
    Particle,
    acceleration,
    kinetic_energy,
    potential_energy,
    total_energy,
)


def main() -> None:
    particle = Particle(
        mass=12.0,
        position=[0.1, 0.0, 0.0],
        velocity=[0.1, 0.0, 0.0],
    )
    potential = HarmonicPotential(
        stiffness=100.0,
        center=[0.0, 0.0, 0.0],
    )

    force = potential.force(particle.position)
    accel = acceleration(particle, force)

    print("Initial state")
    print(f"  mass       = {particle.mass:.3f} Da")
    print(f"  position   = {particle.position} nm")
    print(f"  velocity   = {particle.velocity} nm/ps")
    print()
    print("Energies")
    print(f"  kinetic    = {kinetic_energy(particle):.6f} kJ/mol")
    print(f"  potential  = {potential_energy(particle, potential):.6f} kJ/mol")
    print(f"  total      = {total_energy(particle, potential):.6f} kJ/mol")
    print()
    print("Force and acceleration")
    print(f"  force      = {force} kJ/(mol nm)")
    print(f"  acceleration = {accel} nm/ps^2")
    print()
    print("The negative force points back toward the equilibrium position.")


if __name__ == "__main__":
    main()
