from setuptools import setup
from os import path
#from setuptools_scm import get_version

dir_path = path.abspath(path.dirname(__file__))

with open(path.join(dir_path, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='pdbutil',
    packages=['pdbutil'],
    license='MIT',
    url='https://github.com/ShintaroMinami/pdbutil',
    description='A simple module for handling protein backbone coordinates',
    long_description=long_description,
    long_description_content_type='text/markdown',
    keywords=['pdb', 'backbone'],

    author='Shintaro Minami',
    author_email='shintaro.minami@gmail.com',

    use_scm_version={'local_scheme': 'no-local-version'},

    setup_requires=['setuptools_scm'],
    install_requires=['numpy'],
 
    include_package_data=True,
    scripts=[
    ],

    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.8',
    ],
)