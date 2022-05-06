# RemoveDups
> Removes PCR duplicates from SAM files. 

## Table of Contents
* [General Info](#general-information)
* [Setup](#setup)
* [Usage](#usage)
* [Project Status](#project-status)
* [Acknowledgements](#acknowledgements)



## General Information
- The main purpose of this project is to remove reads which are PCR duplicates from a SAM file. 



## Setup

**Dependencies**

python3.9
> Python packages: argparse, re
matplotlib




## Usage

**Running the script**

```
$ ./remove_dups.py -h
usage: remove_dups.py [-h] -f FILE -u UMI [-p]

Remove PCR duplicates from SAM file and output new SAM file.

optional arguments:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  absolute path to SAM file
  -u UMI, --umi UMI     file containing the list of UMIs (unset if randomers
                        instead of UMIs)
  -p, --paired          flag if file is paired end
```



## Project Status
Project is: _complete_ 


## Acknowledgements
- This project was based on instructions given by the University of Oregon's Knight Campus Internship Program for Bioinformatics and Genomics (BGMP)
- Many thanks to the BGMP instructors and fellow classmates for their feedback



