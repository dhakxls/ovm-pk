"""Example flexible pipeline usage."""
from ovmpk import Pipeline
import yaml

# Load config
with open("configs/flexible_example.yaml") as f:
    config = yaml.safe_load(f)

# Customize for different systems
config["protein"]["pdb"] = "5VCC"  # Try "1A2C", "3UA7", etc.
config["ligand"]["name"] = "ketoconazole"  # Try "fluconazole", "itraconazole"

# Run pipeline
pipeline = Pipeline(config)
results = pipeline.run()
print(f"Pipeline completed. Results: {results}")
