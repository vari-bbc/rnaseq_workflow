#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/.."

runtime="${SINGULARITY_CMD:-}"
if [[ -z "$runtime" ]]; then
    if command -v singularity >/dev/null 2>&1; then
        runtime="singularity"
    else
        echo "Could not find singularity on PATH." >&2
        exit 1
    fi
fi

image_dir="${IMAGE_DIR:-workflow/containers/images}"
r_recipe="${R_CONTAINER_RECIPE:-workflow/containers/renv.def}"
mkdir -p "$image_dir"

targets=("$@")
if [[ "${#targets[@]}" -eq 0 ]]; then
    targets=(rnaseq qc visualisation variants renv)
fi

contains_target() {
    local needle="$1"
    local target
    for target in "${targets[@]}"; do
        [[ "$target" == "$needle" || "$target" == "all" ]] && return 0
    done
    return 1
}

build_image() {
    local name="$1"
    local recipe="$2"
    contains_target "$name" || return 0
    echo "Building $name from $recipe"
    "$runtime" build --force "$image_dir/$name.sif" "$recipe"
}

build_image rnaseq workflow/containers/rnaseq.def
build_image qc workflow/containers/qc.def
build_image visualisation workflow/containers/visualisation.def
build_image variants workflow/containers/variants.def
build_image renv "$r_recipe"
