from setuptools import setup, find_packages

setup(name='SPack_Python',
      version='1.001',
      description='Basic python packages',
      url='https://github.com/hantke/SPack',
      author='Sergio Contreras',
      author_email='el.hantke@gmail.com',
      #license='MIT',
      packages=['SPP_Plot','SPP_Mate','SPP_RW'],
      install_requires=[
          'numpy','matplotlib',
      ],
      zip_safe=False)
#python setup.py build_ext --inplace
#python setup.py install --user
#sudo python setup.py install
