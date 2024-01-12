# Proksee Batch

[![PyPI](https://img.shields.io/pypi/v/proksee-batch.svg)][pypi_]
[![Status](https://img.shields.io/pypi/status/proksee-batch.svg)][status]
[![Python Version](https://img.shields.io/pypi/pyversions/proksee-batch)][python version]
[![License](https://img.shields.io/pypi/l/proksee-batch)][license]

[![Read the documentation at https://proksee-batch.readthedocs.io/](https://img.shields.io/readthedocs/proksee-batch/latest.svg?label=Read%20the%20Docs)][read the docs]
[![Tests](https://github.com/stothard-group/proksee-batch/workflows/Tests/badge.svg)][tests]
[![Codecov](https://codecov.io/gh/stothard-group/proksee-batch/branch/main/graph/badge.svg)][codecov]

[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)][pre-commit]
[![Black](https://img.shields.io/badge/code%20style-black-000000.svg)][black]

[pypi_]: https://pypi.org/project/proksee-batch/
[status]: https://pypi.org/project/proksee-batch/
[python version]: https://pypi.org/project/proksee-batch
[read the docs]: https://proksee-batch.readthedocs.io/
[tests]: https://github.com/stothard-group/proksee-batch/actions?workflow=Tests
[codecov]: https://app.codecov.io/gh/stothard-groupproksee-batch
[pre-commit]: https://github.com/pre-commit/pre-commit
[black]: https://github.com/psf/black

## Features

_Proksee Batch_ is a command-line tool for visualizing multiple prokaryotic
genomes via the [Proksee] web application.

## Installation

You can simply install _Proksee Batch_ via [pip] from [PyPI]:

```console
pip install proksee-batch
```

Or, better yet, use [pipx] to install _Proksee Batch_ from [PyPI] in an isolated environment:

```console
pipx install proksee-batch
```

To install _Proksee Batch_ in a [conda] environment, use a YAML file like this:

```yaml
name: proksee-batch
channels:
  - conda-forge
dependencies:
  - python=3.11
  - pip
  - pip:
      - proksee-batch
```

## Usage

Please see the [Command-line Reference] for details.

## Contributing

Contributions are welcome.
To learn more, see the [Contributor Guide].

## License

Distributed under the terms of the [MIT license][license],
_Proksee Batch_ is free and open source software.

## Issues

If you encounter any problems,
please [file an issue] along with a detailed description.

## Credits

This project was initially generated using [@cjolowicz]'s [Hypermodern Python Cookiecutter] template.

[@cjolowicz]: https://github.com/cjolowicz
[pypi]: https://pypi.org/
[hypermodern python cookiecutter]: https://github.com/cjolowicz/cookiecutter-hypermodern-python
[file an issue]: https://github.com/stothard-group/proksee-batch/issues
[pip]: https://pip.pypa.io/
[pipx]: https://pipx.pypa.io/stable/
[conda]: https://docs.conda.io/en/latest/
[proksee]: https://proksee.ca

<!-- github-only -->

[license]: https://github.com/stothard-group/proksee-batch/blob/main/LICENSE
[contributor guide]: https://github.com/stothard-group/proksee-batch/blob/main/CONTRIBUTING.md
[command-line reference]: https://proksee-batch.readthedocs.io/en/latest/usage.html
