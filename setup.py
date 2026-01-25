from setuptools import setup
from setuptools import find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='CoMBCR',
    version='0.2.1',
    author='Yiping Zou',
    author_email='yipingzou2-c@my.cityu.edu.hk',
    description='A python lib for CoMBCR',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    packages=find_packages(include=['CoMBCR']),
    package_data={'CoMBCR': ['tokenizer/*']},
    include_package_data = True,
    url='https://github.com/deepomicslab/CoMBCR.git/',
    license='MIT',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)

