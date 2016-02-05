from setuptools import setup  #, find_packages

setup(name='pycrispr',
      version='0.1.2',
      description='Crispr Eating and other tools.',
      url='https://github.com/alaindomissy/crispr',
      author='Alain Domissy',
      author_email='alaindomissy@gmail.com',
      license='MIT',
      # packages=find_packages(exclude=["tests"]),
      packages=['crispr', 'pybasespace'],
      install_requires=[
            'biopython==1.66',
            # 'coverage==4.0.3',
            # 'cycler==0.9.0',
            # 'Cython==0.23.4',
            'decorator==4.0.6',
            # 'ipython==4.0.3',
            # 'ipython-genutils==0.1.0',
            # 'nose==1.3.7',

            'numpy==1.10.2',
            'matplotlib==1.5.1',
            # 'pandas==0.17.1',

            'path.py==8.1.2',

            # 'pexpect==3.3',
            # 'pickleshare==0.5',

            'py==1.4.31',
            'pybedtools==0.7.4',

            'pyparsing==2.0.3',
            'pysam==0.8.4',
            # 'pytest==2.8.5',
            'python-dateutil==2.4.2',
            'pytz==2015.7',
            'simplegeneric==0.8.1',
            'six==1.10.0',
            # 'traitlets==4.1.0',
            # 'wheel==0.26.0',
             ],
      classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Utilities",
        "License :: OSI Approved :: MIT License",
      ],
      zip_safe=False)
