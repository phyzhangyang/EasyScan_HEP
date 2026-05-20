from setuptools import setup


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
    packages=["easyscan_hep", "ui"],
    package_data={"ui": ["templates/*.html", "static/*.js", "static/*.css"]},
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
