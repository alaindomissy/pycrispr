from setuptools import setup, find_packages

setup(name='pycrispr',
      version='0.1.2',
      description='Crispr Eating and other tools.',
      url='https://github.com/alaindomissy/buffet',
      author='Alain Domissy',
      author_email='alaindomissy@gmail.com',
      license='MIT',
      # packages=find_packages(exclude=["tests"]),
      packages=['buffet', 'pybasespace'],
      install_requires=[''],
      classifiers=[
        "Development Status :: 4 - Beta",
        "Topic :: Utilities",
        "License :: OSI Approved :: MIT License",
      ],
      zip_safe=False)
