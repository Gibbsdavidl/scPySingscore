from setuptools import setup

def readme():
    with open('../README.rst') as f:
        return f.read()

setup(name='scsingscore',
      version='0.0.1',
      description='modified singscore in Python',
      url='http://github.com/kristyhoran/singscore',
      author='David Gibbs,Kristy Horan',
      author_email='gibbsdavidl@gmail.com, kristyhoran15@gmail.com',
      license='MIT',
      packages=['scsingscore', 'scsingscore.test'],
      install_requires=[
          'pandas', 'numpy', 'matplotlib', 'seaborn', 'scipy', 'statsmodels', 'scanpy', 'tqdm'
      ],
      zip_safe=False)
