from setuptools import setup

setup(
    name='pretty_basemap',
    packages=['pretty_basemap'],
    version='0.0.1',
    description='Avoid Lazy Geopandas Maps',
    author='Sam Koebrich',
    install_requires=[
        'geopandas',
        'cartopy',
        'proj',
        'shapely',
        'pandas',
        'numpy',
        'descartes'
    ]
)
