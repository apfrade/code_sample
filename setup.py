import setuptools

with open("README.md", encoding="utf8") as readme_file:
    long_description = readme_file.read()

requirements = ["pandas>=0.2", "numpy>=1", "pickle"]

setuptools.setup(
    name='ElectrostaticPotential',
    version="0.0.1",
    author='Andre Frade',
    author_email="andre.frade@hertford.ox.ac.uk",
    description='A library of electrostatic potential 2D molecular descriptors',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/apfrade/ConfidenceMeasure.git",
    packages=['electrostatic_potential'],
    packages=setuptools.find_packages(),
    install_requires= requirements,
    python_requires='>=3',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
],  
    # includes non code files 
    packages = find_packages('electrostatic_potential'),
    package_dir = {'':'electrostatic_potential'},
    include_package_data = True,
    # run tests
    test_suite='tests',   
)