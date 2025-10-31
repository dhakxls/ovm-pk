"""Lambda window optimization (Amber-free)."""
from openmm import unit
from openmmtools.alchemy import AbsoluteAlchemicalFactory, AlchemicalRegion

def run_lambda_windows(system, topology, positions, config):
    """Setup lambda windows for alchemical transformations."""
    # Create alchemical region
    alch_region = AlchemicalRegion(
        alchemical_atoms=config["alchemical_atoms"],
        annihilate_electrostatics=True,
        annihilate_sterics=True,
        softcore_alpha=config.get("softcore_alpha", 0.5)
    )
    
    # Create factory
    factory = AbsoluteAlchemicalFactory()
    alchemical_system = factory.create_alchemical_system(system, alch_region)
    
    # Generate lambda schedule
    lambda_schedule = config["lambda_schedule"]
    
    return {
        "systems": [alchemical_system] * len(lambda_schedule),
        "lambdas": lambda_schedule
    }
