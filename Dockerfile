FROM r-base:latest

USER root

RUN apt-get update && apt-get -y upgrade
RUN apt-get install libcurl4-openssl-dev

WORKDIR /tmp

RUN wget https://bioconductor.org/packages/3.3/bioc/src/contrib/RankProd_2.44.0.tar.gz 
RUN wget https://cran.r-project.org/src/contrib/entropy_1.2.1.tar.gz
RUN R -e "install.packages(c('/tmp/RankProd_2.44.0.tar.gz','/tmp/entropy_1.2.1.tar.gz'), repos = NULL, type = 'source')"
RUN R -e "install.packages('factoextra', dependencies = TRUE)"
RUN cd .. && rm -R tmp

RUN mkdir /data
RUN mkdir /home/Entropic_Ranks

ADD Entropic_Ranks_v1.0.R /home/Entropic_Ranks/Entropic_Ranks.R

WORKDIR /home/Entropic_Ranks

#CMD Rscript Entropic_Ranks.R /data/data_table.txt /data/population_vector.txt
