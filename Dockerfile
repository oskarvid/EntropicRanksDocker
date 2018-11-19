FROM r-base:latest

USER root

# RUN apt-get update && apt-get install -y --no-install-recommends apt-utils
# RUN apt-get install apt-utils
# RUN apt-get update && apt-get -y upgrade
# RUN apt-get -f install libcurl4-openssl-dev

# RUN apt-get update && apt-get -y upgrade && apt-get install -y aptitude && aptitude -f install -y libcurl4-openssl-dev

RUN apt-get update && apt-get -y upgrade && apt-get install -y aptitude && aptitude -f install -y libcurl4-openssl-dev && apt-get install -t unstable curl

WORKDIR /tmp

# RUN wget https://bioconductor.org/packages/3.3/bioc/src/contrib/RankProd_2.44.0.tar.gz 
# RUN wget https://cran.r-project.org/src/contrib/entropy_1.2.1.tar.gz
# RUN R -e "update.packages(ask = FALSE)"
# RUN R -e "install.packages(c('curl','ggplot2', 'abind', 'dendextend', 'FactoMineR', 'ggpubr', 'reshape2', 'ggrepel', 'tidyr'), update =F, ask = FALSE)"
# RUN R -e "install.packages(c('/tmp/RankProd_2.44.0.tar.gz','/tmp/entropy_1.2.1.tar.gz'), repos = NULL, type = 'source')"
# RUN R -e "install.packages('factoextra', dependencies = TRUE)"


# RUN R -e "install.packages('libxml-2.0', dependencies = TRUE)"
# RUN R -e "install.packages('xml-2.0', dependencies = TRUE)"
# RUN apt-get install -y libxml2-dev
# RUN aptitude -f install -y libxml2-dev
# RUN apt-get install -y musl-dev
# RUN apt-get install -t unstable curl
# RUN R -e "install.packages('curl', dependencies = TRUE)"
# RUN R -e "install.packages('webutils', dependencies = TRUE)"
# RUN R -e "install.packages('factoextra', dependencies = TRUE)"



RUN cd .. && rm -R tmp

RUN mkdir /data
RUN mkdir /home/Entropic_Ranks
ADD Entropic_Ranks_v1.0.R /home/Entropic_Ranks/Entropic_Ranks.R

WORKDIR /home/Entropic_Ranks

CMD Rscript Entropic_Ranks.R /data/data_table.txt /data/population_vector.txt null 1 FALSE FALSE TRUE TRUE TRUE 2 FALSE
