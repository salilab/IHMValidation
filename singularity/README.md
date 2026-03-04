### Building runtime image

IHMValidation pipeline is designed to work inside the [Apptainer](https://apptainer.org/docs/user/latest/introduction.html) (previously Singularity) container with all 3rd-party dependencies.

To build the image:
1. Navigate to the directory holding `Singularity.def` file:
`cd IHMValidation/singularity`
2. Download [Chimera](https://www.cgl.ucsf.edu/chimera/download.html) and place in the same directory. The current pipeline was tested with `Chimera-1.19` for headless servers.
3. Download [ChimeraX](https://www.cgl.ucsf.edu/chimerax/download.html) `.deb` package for Ubuntu 22.04 and place in the same directory. The current pipeline was tested with `ChimeraX-1.9`.
4. Select an appropriate [timezone](https://en.wikipedia.org/wiki/List_of_tz_database_time_zones). For example `America/Los_Angeles`
5. Build the image:
`apptainer build --build-arg TZ=America/Los_Angeles --build-arg CHIMERA=./chimera-1.19-linux_x86_64_osmesa.bin --build-arg CHIMERAX=./ucsf-chimerax_1.9ubuntu22.04_amd64.deb ihmv.sif Singularity.def`

which will generate `ihmv.cif` image file. 

**NB 1:** To get a date-tagged name, use `ihmv_$(date +%Y%m%d).sif` as an image name. 

**NB 2:** It takes about ~25 minutes on a modern workstation to build the image from scratch. 
