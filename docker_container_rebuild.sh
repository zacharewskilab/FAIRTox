#!/bin/bash

sudo docker stop zed
sudo docker rm zed
sudo docker build -t shinyserver .
sudo docker run -d -v "/home/ec2-user/RData":/srv/shiny-server/RData -p 80:80 --name zed shinyserver
