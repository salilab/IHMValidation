### Building runtime image

IHMValidation pipeline is designed to work inside the [Singularity](https://docs.sylabs.io/guides/4.1/user-guide/introduction.html) container with all 3rd-party dependencies.

To build the image:
1. Navigate to the directory holding `Singularity.def` file:
`cd IHMValidation/docker`
2. Download [ATSAS](https://www.embl-hamburg.de/biosaxs/download.html) `.deb` package for Ubuntu 22.04 and place in the same directory. The current pipeline was tested with `ATSAS-3.2.1` but newer versions should work too.
3. Download [Chimera](https://www.cgl.ucsf.edu/chimera/download.html) and place in the same directory. The current pipeline was tested with `Chimera-1.19` for headless servers.
4. Download [ChimeraX](https://www.cgl.ucsf.edu/chimerax/download.html) `.deb` package for Ubuntu 22.04 and place in the same directory. The current pipeline was tested with `ChimeraX-1.9`.
5. Select an appropriate [timezone](https://en.wikipedia.org/wiki/List_of_tz_database_time_zones). For example `America/Los_Angeles`
6. Build the image:
`singularity build --build-arg TZ=America/Los_Angeles --build-arg ATSAS=./ATSAS-3.2.1-1_amd64.deb --build-arg CHIMERA=./chimera-1.19-linux_x86_64_osmesa.bin --build-arg CHIMERAX=./ucsf-chimerax_1.9ubuntu22.04_amd64.deb ihmv_$(date +%Y%m%d).sif Singularity.def`

which will generate a `.sif` image using current date as a nametag (i.e. `ihmv_20250704.sif`)

**NB 1**: It takes about ~25 minutes on a modern workstation to build the image from scratch. 

**NB 2**: Due to the [Molprobity](https://github.com/rlabduke/MolProbity)'s rolling release model we can't freeze Molprobity version thus we can't guarantee the `IHMValidation` package will work with a freshly-built image. 
