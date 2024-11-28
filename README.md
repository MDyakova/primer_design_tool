# Primer Design Tool for Gene Editing Sites

Welcome to the **Primer Design Tool**, a user-friendly tool for designing primers using BLAST, tailored for gene editing applications. This tool is packaged in a Docker container for ease of deployment and consistent environment setup.

## Features

- **BLAST Integration**: Automates the identification of target sequences for primer design.
- **Amplicon Optimization**: Designs primers with optimal melting temperature (Tm) and GC content for your specific target.
- **Customizable Parameters**: Allows users to specify primer length, Tm range, and product size.
- **Gene Editing Focus**: Targets regions modified by CRISPR, TALENs, or other genome editing technologies.
- **Dockerized**: Simple and consistent deployment across platforms.

## Prerequisites

- Install Docker: [Docker Installation Guide](https://docs.docker.com/get-docker/)

## Getting Started

### 1. Pull the Docker Image
```bash
docker pull mdyakova/primer_design_tool:v1
```

### 2. Run the Docker Container
```bash
docker run -p 5000:5000 -d dockerfile:primerdesign
```

## Usage

See "User_guide.pdf"

## Contact

For questions or suggestions, please contact `m.dyakova.ml@gmail.com`.

---

Happy primer designing! 🎯





## Make new docker image

docker build -t dockerfile:primerdesign ./ 
docker run -p 5000:5000 -d dockerfile:primerdesign

## docker hub
docker tag dockerfile:primerdesign mdyakova/primer_design_tool:v1
docker login
docker push mdyakova/primer_design_tool:v1

## Check code quality
python3 -m black main.py
python3 lint.py --threshold 7
