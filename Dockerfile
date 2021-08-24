FROM python:3
RUN pip3 install pandas
RUN pip3 install numpy
RUN pip3 install networkx
RUN pip3 install argparse

COPY . /app