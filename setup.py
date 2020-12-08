import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="RCAnalysis", # Replace with your own username
    version="0.0.1",
    author="SemiQuant",
    author_email="Jason.Limberis@ucsf.edu",
    description="Takes as input MinION fastq files and primer sequence Searches for the primer sequence in each read cuts them up, creating a fastq file for each.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/SemiQuant/RCAnalysis",
    packages=['Bio.Seq', 'distutils.command', 'Bio', 'fuzzysearch', 'statistics', 'logging', 'argparse'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)
