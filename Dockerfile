FROM grandforest-shiny-base
MAINTAINER Simon J. Larsen <simonhffh@gmail.com>

EXPOSE 3838

RUN sudo apt-get update && sudo apt-get install -y libssl-dev libxml2-dev libcurl4-gnutls-dev

COPY install.R install.R
RUN R -e "source('install.R')"

COPY . /srv/shiny-server

CMD ["/usr/bin/shiny-server.sh"]
