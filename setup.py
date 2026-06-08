from setuptools import setup
from glob import glob
import os


def file_glob(pattern):
    return [path for path in glob(pattern) if os.path.isfile(path)]


CORE_REQUIRES = [
    "numpy>=1.20",
    "scipy>=1.7",
    "matplotlib>=3.5",
    "pandas>=1.3",
]

UI_REQUIRES = [
    "fastapi>=0.95",
    "uvicorn>=0.20",
    "jinja2>=3.0",
    "python-multipart>=0.0.6",
]

DYNESTY_REQUIRES = ["dynesty>=2.1"]
EMCEE_REQUIRES = ["emcee>=3.1"]
MULTINEST_REQUIRES = ["pymultinest>=2.12"]


setup(
    name="easyscan-hep",
    version="2.0",
    description="A tool for connecting programs to scan the parameter space of physics models.",
    long_description=open("README.rst", encoding="utf-8").read(),
    long_description_content_type="text/x-rst",
    python_requires=">=3.9",
    license="Apache-2.0",
    author="Yang Zhang, Liangliang Shang, Yang Xiao",
    url="https://github.com/phyzhangyang/EasyScan_HEP",
    project_urls={
        "Documentation": "https://arxiv.org/pdf/2304.03636.pdf",
        "Source": "https://github.com/phyzhangyang/EasyScan_HEP",
        "Website": "https://easyscanhep.hepforge.org",
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Environment :: Web Environment",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    packages=["easyscan_hep", "ui"],
    package_data={"ui": ["templates/*.html", "static/*.js", "static/*.css"]},
    data_files=[
        ("templates", file_glob("templates/*.ini") + file_glob("templates/*.txt")),
        ("utils", ["utils/TestFunction.py", "utils/TestFunction_input.dat", "utils/OnePointBatch.in"]),
        ("utils/MSSM_mW", [
            "utils/MSSM_mW/for_MSSM_mW.susyhit",
            "utils/MSSM_mW/hbounds.txt",
            "utils/MSSM_mW/makefile",
        ]),
        ("utils/SSM_DM_EWPT", [
            "utils/SSM_DM_EWPT/bootstrap.py",
            "utils/SSM_DM_EWPT/plot_dm_ewpt.py",
            "utils/SSM_DM_EWPT/README.md",
        ]),
        ("utils/SSM_DM_EWPT/patches", [
            "utils/SSM_DM_EWPT/patches/micromegas_7.1_easyscan.patch",
        ]),
    ],
    install_requires=CORE_REQUIRES,
    extras_require={
        "ui": UI_REQUIRES,
        "dynesty": DYNESTY_REQUIRES,
        "emcee": EMCEE_REQUIRES,
        "multinest": MULTINEST_REQUIRES,
        "samplers": DYNESTY_REQUIRES + EMCEE_REQUIRES,
        "all": UI_REQUIRES + DYNESTY_REQUIRES + EMCEE_REQUIRES,
    },
    py_modules=[
        "auxfun",
        "config_checker",
        "constraint",
        "initialize",
        "ploter",
        "program",
        "readin_config",
        "scan_controller",
        "scanner",
        "statfun",
    ],
    package_dir={"": "src", "easyscan_hep": "easyscan_hep", "ui": "ui"},
    entry_points={"console_scripts": ["easyscan=easyscan_hep.cli:main"]},
)
