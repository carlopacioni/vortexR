FROM ghcr.io/rocker-org/devcontainer/r-ver:4.2

# To use a more complex devcontainer, modify devcontainer.json to use this Dockerfile 
# by replacing the line "image" with "dockerFile": "Dockerfile"

# non interactive frontend for locales
ENV DEBIAN_FRONTEND=noninteractive

# Install system packages
# RUN apt-get update && \
#    apt-get -y install --no-install-recommends \
#    PACKAGES HERE  && \
#    rm -rf /var/lib/apt/lists/*

# generating locales
# RUN sed -i -e 's/# en_GB.UTF-8 UTF-8/en_GB.UTF-8 UTF-8/' /etc/locale.gen && \
#     dpkg-reconfigure --frontend=noninteractive locales && \
#     update-locale LANG=en_GB.UTF-8
# ENV LANGUAGE=en_GB.UTF-8 LANG=en_GB.UTF-8 LC_ALL=en_GB.UTF-8

# installing cpanm & missing latexindent dependencies
# RUN curl -L http://cpanmin.us | perl - --self-upgrade && \
#     cpanm Log::Dispatch::File YAML::Tiny File::HomeDir
