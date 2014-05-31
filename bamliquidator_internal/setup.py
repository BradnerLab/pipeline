from setuptools import setup

setup(
    name='BamLiquidatorBatch',
    version='0.6.14',
    maintainer='John DiMatteo',
    maintainer_email='jdimatteo@gmail.com',
    packages=['bamliquidatorbatch'],
    url='https://github.com/BradnerLab/pipeline/wiki/bamliquidator',
    license='LICENSE.txt',
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
