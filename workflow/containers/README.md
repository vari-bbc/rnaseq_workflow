# Singularity images

The workflow can run with local SIF images or hosted image URIs. The default
configuration points to local SIF paths under `workflow/containers/images/`.
Build them from the repository root with:

```bash
bin/build_singularity_images.sh
```

The R image is intentionally not conda based. `renv.def` starts from Rocker R,
installs the system libraries used by the workflow, and preinstalls the packages
listed in `config/R_proj_packages.txt` plus common organism databases into the
image site library. This image can be large; after building it, hosting it in a
registry or shared image store is usually better than rebuilding it per project.

For a smaller image that lets `renv` restore packages at runtime, build with:

```bash
R_CONTAINER_RECIPE=workflow/containers/renv-minimal.def bin/build_singularity_images.sh renv
```

To use hosted images instead of local SIFs, replace the entries under
`containers:` in `config/config.yaml` with URIs such as `docker://...`,
`oras://...`, or a shared absolute `.sif` path.
