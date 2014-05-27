from distutils.core import setup

setup(
    name='BamLiquidator',
    packages=['bamliquidator'],
    install_requires=[
        "tables >= 3.0.0",
        "bokeh >= 4.0.0",
    ],
)
