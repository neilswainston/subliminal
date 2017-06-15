FROM python:2.7

# Make current directory visible inside Docker container:
COPY . /
WORKDIR /

# Install / update relevant ubuntu packages:
RUN apt-get update \
	&& apt-get install -y --no-install-recommends build-essential curl libgmp3-dev python-pip

# Download and install glpk:
RUN mkdir /usr/local/glpk \
	&& curl http://ftp.gnu.org/gnu/glpk/glpk-4.39.tar.gz \
	| tar xvzC /usr/local/glpk --strip-components=1 \
	&& cd /usr/local/glpk \
	&& ./configure \
	&& make \
	&& make install
	
# Install pyglpk:
RUN pip install --upgrade pip \
	&& pip install setuptools \
	&& pip install glpk

# Updata paths:
ENV LD_LIBRARY_PATH /usr/local/lib:${LD_LIBRARY_PATH}

# Run test:
ENTRYPOINT ["python"]
CMD ["test.py"]