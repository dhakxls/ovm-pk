from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="ovm-pk",
    version="0.1.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="Physics-Aware Docking & Validation Pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/ovm-pk",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.8",
    install_requires=[
        "click>=8.0.0",
        "requests>=2.25.0",
        "numpy>=1.20.0",
        "pandas>=1.2.0",
        "matplotlib>=3.3.0",
        "rich>=10.0.0",
        "questionary>=1.10.0",
        "openmm>=7.6.0",
        "mdtraj>=1.9.0",
        "rdkit>=2021.03.0",
    ],
    entry_points={
        "console_scripts": [
            "ovmpk=ovmpk.cli.main:cli",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
)
