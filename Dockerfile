# ============================================================
# Dockerfile — builds a container image that runs forest
# simulation API + front-end as a single, self-contained app
#
# local environment:
#   R version 4.4.0 (aarch64-apple-darwin20 / Apple Silicon)
# ============================================================

# ---- 1. FROM ----------------------------------
#
# NOTE: pinned to 4.5.2 rather than your exact local 4.4.0 — Docker Hub
# confirms 4.4.0 only has an amd64 build published, no arm64 variant
FROM rocker/tidyverse:4.5.2

# ---- 2. RUN: install system libraries ----------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    libsodium-dev \
    && rm -rf /var/lib/apt/lists/*

# ---- 3. RUN: install the R packages plumber.R needs -------
# tidyverse (dplyr, purrr, tibble, etc.) is already in the base image
# need to add:
#   plumber - runs the API itself
#   RANN    - k-d tree nearest-neighbor search, used in your shading calc
RUN R -e "install.packages(c('plumber', 'RANN'), repos = 'https://cloud.r-project.org')"

# ---- 4. WORKDIR: set the working directory inside the container --------
WORKDIR /app

# ---- 5. COPY: bring your project files into the image -------------------
# Order matters here 
COPY R/ ./R/
COPY data/ ./data/
COPY www/ ./www/
COPY plumber.R ./

# ---- 6. EXPOSE: document which port the app listens on -------------------
EXPOSE 8000

# ---- 7. CMD: what to run when a container starts -------------------------
CMD ["R", "-e", "plumber::plumb('plumber.R') |> plumber::pr_run(host='0.0.0.0', port=8000)"]
