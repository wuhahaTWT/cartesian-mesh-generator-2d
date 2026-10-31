FROM gas-turbine-hydrogym:1.0.0-fd2026.4.2
# Existing Gmsh Python wheel only lacks its system OpenGL utility dependency.
RUN apt-get update && apt-get install -y --no-install-recommends libglu1-mesa libgl1 libxft2 \
    && rm -rf /var/lib/apt/lists/*
ENV OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg PYTHONUNBUFFERED=1
