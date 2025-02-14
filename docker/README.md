# Building runtime image

We switched from Docker to [Singularity](https://docs.sylabs.io/guides/4.1/user-guide/introduction.html) containers to run the validation pipeline. 

To build the image:
1. Navigate to the directory holding `Singularity.def` file:
`cd IHMValidation/docker`
2. Download [ATSAS](https://www.embl-hamburg.de/biosaxs/download.html) `.deb` package for Ubuntu 22.04 and put into the same directory. The current pipeline was tested with `ATSAS-3.2.1` but newer versions should work too.
3. Select an appropriate [timezone](https://en.wikipedia.org/wiki/List_of_tz_database_time_zones). For example `America/Los_Angeles`
4. Build the image:
`singularity build --build-arg TZ=America/Los_Angeles --build-arg ATSAS=./ATSAS-3.2.1-1_amd64.deb ihmv_$(date +%Y%m%d).sif Singularity.def`
which will generate a `.sif` image using current date as a nametag (i.e. `ihmv_20240227.sif`)

**NB 1**: It takes about ~25 minutes on a modern workstation to build the image from scratch. 

**NB 2**: Due to the [Molprobity](https://github.com/rlabduke/MolProbity)'s rolling release model we can't freeze Molprobity version thus we can't guarantee the `IHMValidation` package will work with a freshly-built image. 
