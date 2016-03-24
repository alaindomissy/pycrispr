from setuptools import setup  #, find_packages

setup(name='crispr',
      version='0.1.2',
      description='Crispr Eating and other tools.',
      url='https://github.com/alaindomissy/crispr',
      author='Alain Domissy',
      author_email='alaindomissy@gmail.com',
      license='MIT',
      # packages=find_packages(exclude=["tests"]),
      packages=['crispr'],
      install_requires=[

            #'numpy==1.10.4',      # gets installed by conda, would need gcc if done here
                                   # 1.10.2 would be a downgrade from what conda has: 1.10.4
            # 'pybedtools==0.7',     # gets installed by conda

            # 'cycler==0.9.0',
            # 'Cython==0.23.4',
            # 'decorator==4.0.6',

            'biopython==1.66',

            'matplotlib==1.5.1',
            'path.py==8.1.2',

            'pyparsing==2.0.3',
            'pysam==0.8.4',

            'python-dateutil==2.4.2',
            'pytz==2015.7',

            'six==1.10.0',

            'basespaceapp',

            'primer3-py'       # needed for travis CI

            ],
      classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Utilities",
        "License :: OSI Approved :: MIT License",
      ],
      zip_safe=False)
