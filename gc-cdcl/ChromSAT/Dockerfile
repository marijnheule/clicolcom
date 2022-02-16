FROM centos:7

RUN yum -y install yum install centos-release-scl && \
    yum-config-manager --enable rhel-server-rhscl-7-rpms && \
    yum -y install make git ninja mercurial \
    	   centos-release-scl devtoolset-7 rh-python36 python-pip \
	   centos-release-scl zlib-devel boost boost-devel && \
    scl enable rh-python36 "pip install setuptools wheel  && pip install meson"

RUN yum -y install \
    https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm && \
    yum -y install ninja-build

WORKDIR /tmp/build/

RUN hg clone http://hg@bitbucket.org/gkatsi/minicsp

#ADD ./tclap /tmp/build/minicsp/tclap
#
ADD ./src /tmp/build/minicsp/src
ADD ./tclap /tmp/build/minicsp/tclap
ADD ./sparsehash /tmp/build/minicsp/sparsehash
ADD ./sota /tmp/build/minicsp/sota
ADD ./meson.build /tmp/build/minicsp
ADD ./meson_options.txt /tmp/build/minicsp

RUN scl enable devtoolset-7 rh-python36 "meson --buildtype=debugoptimized minicsp minicsp/debugopt" && scl enable devtoolset-7 "ninja-build -C minicsp/debugopt/"

RUN ln -s /tmp/build/minicsp/debugopt/gc /usr/bin
