from setuptools import setup

setup(
    name='speakeasy',
    version='0.1',
    description='The speakeasy tool automates conversion of SMIRNOFF parameters into AMBER format.',
    url='https://github.com/nhenriksen/speakeasy',
    author='Niel M. Henriksen',
    author_email='shireham@gmail.com',
    license='MIT',
    packages=['speakeasy'],
    zip_safe=False,

    # Create command line script
    entry_points={
        'console_scripts': [
            'speakeasy = speakeasy.command:main',
        ],
    },

)

