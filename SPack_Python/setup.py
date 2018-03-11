from setuptools import setup, find_packages

setup(name='SPack_Python',
      version='1.001',
      description='Basic python packages',
      url='https://github.com/hantke/SPack',
      author='Sergio Contreras',
      author_email='el.hantke@gmail.com',
      #packages=find_packages(),
      package_dir = {
            'SPack_Python': 'SPack_Python',
            'SPack_Python.Plot': 'SPack_Python/Plot',
            'SPack_Python.Mate': 'SPack_Python/Mate',
            'SPack_Python.RW': 'SPack_Python/RW',
            'SPack_Python.Advance': 'SPack_Python/Advance'},
      packages=['SPack_Python', 'SPack_Python.Plot', 'SPack_Python.Mate',
                  'SPack_Python.RW','SPack_Python.Advance'],
      install_requires=[
          'numpy','matplotlib'
      ],
      zip_safe=False)
#python setup.py build_ext --inplace
#python setup.py install --user
