from setuptools import setup

setup(
    name='BamLiquidatorBatch',
    version='0.6.15',
    maintainer='John DiMatteo',
    maintainer_email='jdimatteo@gmail.com',
    packages=['bamliquidatorbatch'],
    url='https://github.com/BradnerLab/pipeline/wiki/bamliquidator',
    license='The MIT License (MIT)',
    entry_points = {
        'console_scripts': [
            'bamliquidator_batch = bamliquidatorbatch.bamliquidator_batch:main',
            'bamliquidator_flattener = bamliquidatorbatch.flattener:main'
        ]
    },
    install_requires=[
        "bokeh >= 0.4.0"
    ],
)
