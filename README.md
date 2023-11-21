# Rowan Python Library

[![pypi](https://img.shields.io/pypi/v/rowan-python.svg)](https://pypi.python.org/pypi/rowan-python)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/charliermarsh/ruff/main/assets/badge/v1.json)](https://github.com/charliermarsh/ruff)

The Rowan Python library provides convenient access to the Rowan API from applications written in the Python language.

## Documentation

The documentation is available [here](https://docs.rowansci.com).

## Installation

To install, run `pip install rowan-python`.

## Usage

Rowan can be run in either blocking (wait until job is complete) or non-blocking (don't wait) modes.
Both modes require generation of an API key at [labs.rowansci.com](https://labs.rowansci.com/account/api-keys).

For now, molecules are specified through [*cctk*](https://cctk.rtfd.io). Additional ways to specify molecules will be added in the future.

Results are returned as dictionaries in the [*stjames*](https://github.com/rowansci/stjames) format.

### Blocking
```
import cctk
import rowan

rowan.api_key = "rowan-sk..."
client = rowan.Client()

# load molecule by name (cctk can also load in a variety of file formats)
molecule = cctk.Molecule.new_from_name("cyclobutane")

# run calculation remotely and return result
result = client.compute(molecule, name="opt cyclobutane", method="b97-3c", tasks=["optimize", "charge"])
print(result)
```

### Non-Blocking
```
import cctk
import rowan

rowan.api_key = "rowan-sk..."
client = rowan.Client(blocking=False)

# load molecule by name (cctk can also load in a variety of file formats)
molecule = cctk.Molecule.new_from_name("cyclobutane")

# start calculation and return id
calc_id = client.compute(molecule, name="opt cyclobutane", method="b97-3c", tasks=["optimize", "charge"])

# retrieve result (and status)
result = client.get(calc_id)
print(result)

# alternately, cancel queued or running calculation
client.stop(calc_id)
```

## Issues

To report issues, please use the "Issues" tab above.

*Corin Wagen, 2023*
