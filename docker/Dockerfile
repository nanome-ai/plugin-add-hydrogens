FROM nanome/nanobabel
ENV ARGS=""

RUN apt install software-properties-common
RUN add-apt-repository -y ppa:deadsnakes/ppa
RUN apt update
RUN apt install -y python3.7
RUN apt install -y python3-pip
RUN apt install -y git

RUN pip3 install --upgrade pip
RUN git clone https://github.com/nanome-ai/nanome.git
RUN cd nanome && git checkout BetterPluginIntegration && python3.7 -m pip install .
# RUN python3.7 -m pip install nanome

COPY . /app
WORKDIR /app

CMD python3.7 run.py ${ARGS}