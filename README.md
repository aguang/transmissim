[![Build Status](https://travis-ci.org/aguang/transmissim.svg?branch=master)](https://travis-ci.org/aguang/transmissim)

# Overview

transmissim is a pipeline that simulates transmission networks, transmission trees, viral genomes evolving along transmission trees, and high throughput sequencing reads. It uses outbreaker, pyvolve, and ART at the moment.

# Usage

To run:

	python simulate.py -p params.txt

Different options for simulation parameters can be edited in params.txt.

# Pipeline Organization

## Transmission Network

Transmission network simulation is done using [outbreaker](https://sites.google.com/site/therepiproject/r-pac/outbreaker).

## Transmission Tree

The transmission tree is currently a binary representation of the transmission network with a python script. We will incorporate coalescent events on the tree in the future.

## Sequence Simulation

Sequence simulation along the transmission tree is done using [pyvolve](https://github.com/sjspielman/pyvolve).

## Read Simulation

Read simulation for each sequence, sorted by taxa is done using [ART](http://www.niehs.nih.gov/research/resources/software/biostatistics/art/).
