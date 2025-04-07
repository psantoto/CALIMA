# CALIMA: Component Analysis and Line Identification for MAsers
CALIMA (Component Analysis and Line Identification for MAsers) is a Python script that uses
[CASA](https://casa.nrao.edu/) functionalities to detect, isolate and quantify
different spectral components in interferometric data cubes of maser
observations.

The latest version is v1.0 (2025).


![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)
![Version](https://img.shields.io/github/v/release/psantoto/CALIMA)
![Last commit](https://img.shields.io/github/last-commit/psantoto/CALIMA)
![Open Issues](https://img.shields.io/github/issues/psantoto/CALIMA)
![Pull Requests](https://img.shields.io/github/issues-pr/psantoto/CALIMA)

---

## Author
Pablo Santo-Tomas

Instituto de Astrofísica de Andalucía (IAA-CSIC)

ORCID: [0009-0003-5228-2804](https://orcid.org/0009-0003-5228-2804)

Email: psantoto@iaa.es

## Description
CALIMA uses CASA's built-in functions to identify the spectral components of a
source in an interferometric image. Although it was developed for the detection
of components of **maser** emission, it can be used on any unresolved
(point-like) line emission. Given a FITS file with a spectral dimension
(i.e. a data cube), CALIMA will return the spectral channels at which the
source's emission peaks, as well as the spatial coordinates of these components.
Currently, CALIMA only works if the two first axes of the data cube are spatial,
while the frequency axis is either the third or fourth one.

Note, however, that CALIMA is not a source-finding algorithm: the spatial region
within the image in which to search for components must be specified. CALIMA then
looks for spectral components in that region and performs a Gaussian fit at each
spectral channel where a component has been found, which returns its spatial
coordinates and the corresponding errors. A more in-depth description on how the
script functions is found in [How it works](#How-it-works).

## Installation

### Requirements
Before running CALIMA, make sure you have:
+ [CASA installed](https://casa.nrao.edu/casa_obtaining.shtml) (only versions 6.6
and 6.7 tested) and properly set up.
+ The [pandas](https://pandas.pydata.org/) library installed in CASA's Python
environment. If this is not the case, CALIMA will install it using pip.

### Download CALIMA
Clone the repository:

```bash
git clone https://github.com/psantoto/CALIMA.git
```

Or download the ZIP:

```bash
wget https://github.com/psantoto/CALIMA/archive/refs/heads/main.zip -O CALIMA.zip
unzip CALIMA.zip
```

## Usage
To run CALIMA, simply start CASA and execute the script:

```python
execfile('path/to/CALIMA.py')
```

You will then be asked to write the name of the file from which CALIMA will take
the parameters. If this file is in a different directory, write the full path to
it. An example parameter file is included in the repository
([parameter-file.txt](./parameter-file.txt)). All parameters are explained in
[How it works](#How-it-works).

Once the script has finished running, an output .csv file will be saved in your
current directory. It includes information about the detected components
(position, frequency/velocity, flux), as well as a header with information about
the image (reference pixel values, increments). An example output file has been
included in the repository ([output-file.csv](./output-file.csv)).

## How it works

### Input
The first thing that CALIMA does is reading the input file given by the user.
The format of this file is rather simple: each parameter must be given in
a different line, each line having the key (name) assigned to the parameter
and its value separated by a colon
(see [parameter-file.txt](./parameter-file.txt)). Most input variables have a
default value, which is used in case it is not specified by the user, or in
case the program is not able to read the input file.

### Detection
Once the parameters have been read, CALIMA operates in two steps: the detection,
which results in a list of spectral channels with emission, and the fitting,
which returns the spatial coordinates of the components previously found.
For the detection, CALIMA first creates a list with the maximum intensity at each
spectral channel. To make sure this maxima correspond to the source, they are
not searched for in the whole image (given by `File Name`), but in a box region
specified by the user: `Source box`. For the program to work properly, this box
must contain the (point-like) source in all channels, while being small enough
to avoid as much noise as possible. After creating this list, a similar one is
generated with the RMS (Root Mean Square) at each channel in a box region
specified by `Noise box`. In this case, the box where the noise is calculated
must not include the source and it should be big enough so that the RMS in it
serves as a good representation of the noise level in the whole image (box area
larger than approximately 10 times the beam solid angle). Both box parameters
must be written in the input file following CASA's format for box regions: if
`x` and `y` are the bottom left corner coordinates and `X`and `Y` are the top
right corner coordinates of the box, then the whole region is specified as
`Source/Noise box: x, y, X, Y`. Note that these coordinates mus be written
in pixel units.

With these two lists, CALIMA then runs a cycle to check for the first of the
three conditions that a maximum must fulfill in order to be considered a
component. This is, an element of the list is considered a candidate for
component when its signal-to-noise ratio is at least `Nsigma source` (default is
3) and it is a local maximum among the maxima in the list. In noisy data, local
maxima can be found everywhere, since the noise creates fluctuations in the
spectrum. To account for this effect, a maximum of the list is considered a
local maximum when it is either the highest valued element of the list (total
maximum) or there is a valley before the next higher valued element is reached.
A valley is an element of the list that is at a distance (in flux) of
`Nsigma source` sigmas below the previous local maxima. Therefore, candidates
are always separated by valleys. If there is not a valley between two apparent
local maxima, only the higher valued one will be considered as such. This
definition implies that all elements valued zero or bellow will always be
valleys. However, the noise in the data induces an error on what we can consider
as zero emission. To account for that uncertainty, CALIMA considers that any
element valued below `Nsigma valley` sigmas (default is 2) is equivalent to
zero, and therefore it is automatically a valley. This condition is especially
important for the detection of components whose signal-to-noise is close to
`Nsigma source`. It is mandatory that `Nsigma valley` be smaller than
`Nsigma source`.

Once this cycle has ended, the resulting candidates must fufill a second
condition in order to be considered components. This is, candidates must have a
total of `N adjacent chans` adjacent channels with emission above
`Nsigma adjacent` sigmas (default for the first parameter is 2, default for the
latter is 3). These count of `N adjacent chans` does not include the channel
with the candidate spectral peak. We consider consecutive channels either
before, around or after the candidate. This second condition is essential in
order to avoid false detections.

As an example, if the candidate is in the spectral channel 60 and
`N adjacent chans` is set to 2 with `Nsigma adjacent: 3`, then there are three
posibilities for the candidate to fulfill this condition:
1. maxima in channels 58 and 59 are valued over 3 sigma,
2. maxima in channels 59 and 61 are valued over 3 sigma or
3. maxima in channels 61 and 62 are valued over 3 sigma.

If we include the candidate, in all three cases there will be 3 (2+1)
consecutive channels with signal over 3 sigma (since `Nsigma source` should be
set at a value equal or higher than `Nsigma adjacent`).

The resulting candidates have to fulfill one final condition before being
considered components. It can often happen that the consecutive maxima found in
the previous condition do not correspond to the same source (i.e. some of the
maxima may actually be produced by noise), even considering the fact that we are
restrained to a box around the source. To make sure these adjacent maxima
correspond to the same source as the candidate, CALIMA places an ellipse shaped
as the beam and centered around the candidate. It then checks whether each
considered adjacent maxima was detected inside this ellipse or not. If not, it
is discarded, which can lead to the candidate failing the previous condition and
being dropped. The size of the ellipse can be changed with the parameter
`Ellipse size`. The default is 0.33, this is, each axis of the ellipse is one
third of the size of the corresponding full width at half maximum of the beam.

### Fitting
The candidates fulfilling all three conditions are considered components. CALIMA
then performs a Gaussian fit at each channel where a component has been found.
The fit is restricted to `Source box` and uses
[CASA's imfit task](https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.analysis.imfit.html).
As a result of these fits, the spatial coordinates of each maximum are returned,
along with the corresponding errors. For the current version of CALIMA, only
point-like sources with a single (spatial) component are supported for this
fitting.

### Output
CALIMA returns all information in a single output file. This file contains a
header with the reference pixel of the image, its reference values for all axes
(position, flux and frequency/velocity) and the increment per pixel in
each of these axes. Then a table in csv format (comma-separated values) shows
the details of each detected component:
+ the channel,
+ the corresponding frequency and velocity,
+ the intensity of the peak and its corresponding error,
+ the flux density over the fitted component and its corresponding error,
+ the position (both in pixels and in angle units) with the corresponding errors.

The name of the ouput file is specified by `Output file` in the parameter file.
An extra file with more information about the fitting of each component can be
created if `Summary flag ` is set to True (default is False). Moreover, CALIMA
does not erase any variables used during the run, meaning that, as long as the
current CASA session is not finished, the user has access to all of them.
Finally, if `Verbose` is set to True (which is the default value), CALIMA will
print on the Terminal the names of all output files, as well as the channels
where the components were detected.

## Contributing
Thank you for considering contributing to CALIMA! If you encounter any issues,
have suggestions, or want to contribute code, please feel free to reach out
using one of the following methods:

1. **Create an Issue on GitHub**
     If you find a bug, have a feature request, or need help, please
[open an issue](https://github.com/psantoto/CALIMA/issues/new) in the
GitHub repository. This allows others to see the issue and contribute to the
discussion.

2. **Submit a Pull Request on GitHub**
     Feel free to
[create a pull request](https://github.com/psantoto/CALIMA/pulls) if you
want to improve the code.

3. **Contact via Email**
     For direct inquiries, bug reports, or suggestions, you can reach
me at: [psantoto@iaa.es](mailto:psantoto@iaa.es).

## How to cite CALIMA
If you use CALIMA in your research, please **cite it properly** to acknowledge
this work. You can find the citation information in [CITATION.cff](./CITATION.cff).

If you use CALIMA in your reasearch and publish results, you **must cite the
software** as follows:

**APA**
> Santo-Tomas, P. (2025). *CALIMA: Component Analysis and Line Identification for MAsers* (Version 1.0) [Software]. Zenodo. DOI: 10.5281/zenodo.15166803

**BibTex**
```bibtex
@software{calima,
    author = {Pablo Santo-Tomas},
    title = {CALIMA: Component Analysis and Line Identification for MAsers},
    year = {2025},
    version = {1.0},
    doi = {10.5281/zenodo.15166803},
    url = {https://github.com/psantoto/CALIMA}
}
```

## Credits
### :uk: English
The development of this code was supported by project ref. AST22_00001_Subp X with funding from the European Union - NextGenerationEU»; the Spanish Ministry of Science, Innovation and Universities; the Spanish Recovery, Transformation and Resilience Plan; the Department of University, Research and Innovation of the Andalusian Regional Government and Consejo Superior de Investigaciones Científicas. We also acknowledge support from grants PID2020-114461GB-I00, PID2023-146295NB-I00, and CEX2021-001131-S, funded by MCIU/AEI/10.13039/501100011033.

### :es: Español
Financiado por el proyecto ref. AST22_00001_Subp X con financiación de la Unión Europea - NextGenerationEU»; el Ministerio de Ciencia, Innovación y Universidades ; el Plan de Recuperación, Transformación y Resiliencia ; la Consejería de Universidad, Investigación e Innovación de la Junta de Andalucía y el Consejo Superior de Investigaciones Científicas. También reconocemos el apoyo de las ayudas PID2020-114461GB-I00, PID2023-146295NB-I00 y CEX2021-001131-S, financiadas por MCIU/AEI/10.13039/501100011033.


<p align="center">
    <img src="images/EN_Funded_by_the_European_Union_POS.jpg" alt="EU logo" style="width:18%;">
    <img src="images/MCIU.Gob.Web.jpg" alt="MCIU logo" style="width:19%;">
    <img src="images/Logo_PRTR_vertical_COLOR.jpg" alt="PRTR logo" style="width:18%;">
    <img src="images/LogoJunta.png" alt="Junta de Andalucia logo" style="width:15%;">
    <img src="images/CSIC.jpg" alt="CSIC logo" style="width:20%;">
</p>

## License
CALIMA is licensed under the **GNU Lesser General Public Licence (GNU GPL v3)**,
which **requires attribution**.  You are free to share and adapt the material as
long as you give appropriate credit and distribute your contributions under the
same license.  If you distribute or modify CALIMA, you must:
+ **Credit the original authors**
+ **Provide a copy of the license**  and the **original copyright** in your
distribution (unless used as a library)
+ **Indicate any modifications**


For more details, see the [LICENSE](LICENSE) file or visit the [official license page](https://www.gnu.org/licenses/lgpl-3.0.html).
