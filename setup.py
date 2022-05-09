from setuptools import setup, find_packages

setup(
    name="genemethods",
    version="0.0.0.53",
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'SequenceExtractor = genemethods.SequenceExtractor.src.sequenceExtractor:main',
        ],
    },
    include_package_data=True,
    author="Adam Koziol",
    author_email="adam.koziol@inspection.gc.ca",
    url="https://github.com/OLC-Bioinformatics/genemethods",
)
