ARG BASE_IMAGE_TAG=latest
FROM quay.io/qiime2/tiny:${BASE_IMAGE_TAG}

ARG VERSION

COPY . /app

WORKDIR /app

# replace "conda-env" within the environment.yml with qiime2-tiny-$VERSION
RUN sed -i "s/conda-env/qiime2-tiny-${VERSION}/g" /app/environment.yml

RUN conda install -n base -c conda-forge mamba

# Install dependencies using mamba and pip
RUN mamba env update -n qiime2-tiny-${VERSION} -f environment.yml

# Install the plugin
RUN mamba run -n qiime2-tiny-${VERSION} pip install .

# Refresh QIIME cache
RUN mamba run -n qiime2-tiny-${VERSION} qiime dev refresh-cache

WORKDIR /data

# Clean up
RUN rm -r /app
