from setuptools import setup

# todo: specify the license in setup

setup(
    name='BamLiquidatorBatch',
    version='0.6.13',
    maintainer='John DiMatteo',
    maintainer_email='jdimatteo@gmail.com',
    packages=['bamliquidatorbatch'],
    url='https://github.com/BradnerLab/pipeline/wiki/bamliquidator',
    entry_points = {
        'console_scripts': [
            'bamliquidator_batch = bamliquidatorbatch.bamliquidator_batch:main'
        ]
    },
    install_requires=[
        "bokeh >= 0.4.0"
    ],
)
