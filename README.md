# sagepy
A python interface to the [SAGE](https://github.com/lazear/sage) search engine for mass spectrometry proteomics.

This repository hosts the main codebase for the sagepy project, which is dedicated to creating a fully functional Python interface for the powerful Sage search engine, originally written in Rust.

The project is structured as follows:

* `sagepy-connector`: This crate creates a Python interface using [PyO3](https://github.com/PyO3) to bind Rust to Python.
* `sagepy`: A pure Python, fully Pythonic wrapper around the exposed Rust code.
*	`qfdrust`: This crate implements basic false discovery rate (FDR) estimation using TDC, following the methods proposed by [Crema](https://github.com/Noble-Lab/crema).
*	`unimod`: A work-in-progress crate that bridges Sage-style PSM annotation with the UNIMOD standard.

There is still much work to be done, but the essentials for peptide identification and FDR estimation are already fully functional.

Get started quickly by installing sagepy via pip:
```
pip install sagepy
```

Do you have any questions or want to contribute? Feel free to reach out at any time!
