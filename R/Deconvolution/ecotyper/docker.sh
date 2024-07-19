#!/usr/bin/env bash
#' 1. Setup docker on terminal and add non-root user:
sudo snap install docker
sudo groupadd docker
sudo usermod -aG docker $researcher
newgrp docker
docker run hello-world
docker pull cibersortx/fractions
docker pull cibersortx/hires
