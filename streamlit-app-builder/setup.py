from setuptools import setup, find_packages

setup(
    name='streamlit-app-builder',
    version='0.1',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    include_package_data=True,
    install_requires=[
        'streamlit',
        'pandas',
        'numpy',
        'matplotlib',
        'plotly',
        'requests',
        'pyteomics',
        'xlsxwriter',
        'spectrum_utils',
    ],
    entry_points={
        'console_scripts': [
            'ms2_main=ms2_main:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: Microsoft :: Windows',
    ],
)