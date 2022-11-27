from setuptools import find_packages, setup  # pragma: no cover

setup(  # pragma: no cover
    name="cmc-broswer",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    author="Peter Smith",
    email="smith.peter.902@gmail.com",
    license="MIT",
    description="Interact with CMC models on disk",
    copyright="Copyright 2022 Peter Smith",
    python_requires=">=3.10",
    version="0.1.0",
    include_package_data=True,
)
