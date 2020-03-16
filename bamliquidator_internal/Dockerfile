FROM ubuntu:bionic as builder

RUN apt-get -y update && DEBIAN_FRONTEND=noninteractive apt-get install -y libbam-dev libhdf5-serial-dev libboost-dev \
    libboost-timer-dev libgoogle-perftools-dev libtbb-dev samtools build-essential

COPY . /opt/liquidator

WORKDIR /opt/liquidator
RUN make -j6


FROM ubuntu:bionic as runner

RUN apt-get -y update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends python3 \
    python3-tables python3-scipy  \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/liquidator
COPY --from=builder /usr/lib/x86_64-linux-gnu/libtcmalloc_minimal.so.4 \
                    /usr/lib/x86_64-linux-gnu/libtbb.so.2 \
                    /usr/lib/x86_64-linux-gnu/libsz.so.2 \
                    /usr/lib/x86_64-linux-gnu/libhdf5_serial.so.100 \
                    /usr/lib/x86_64-linux-gnu/libhdf5_serial_hl.so.100 \
                    /usr/lib/x86_64-linux-gnu/libboost_system.so.1.65.1 \
                    /usr/lib/x86_64-linux-gnu/libaec.so.0 \
                    /usr/lib/x86_64-linux-gnu/libhts.so.2 \
                    /usr/lib/x86_64-linux-gnu/libcurl-gnutls.so.4 \
                    /usr/lib/x86_64-linux-gnu/libnghttp2.so.14 \
                    /usr/lib/x86_64-linux-gnu/librtmp.so.1 \
                    /usr/lib/x86_64-linux-gnu/libpsl.so.5 \
                    /usr/lib/x86_64-linux-gnu/libgssapi_krb5.so.2 \
                    /usr/lib/x86_64-linux-gnu/libkrb5.so.3 \
                    /usr/lib/x86_64-linux-gnu/libk5crypto.so.3 \
                    /usr/lib/x86_64-linux-gnu/libkrb5support.so.0 \
                    /lib/x86_64-linux-gnu/libkeyutils.so.1 \
                    /usr/lib/x86_64-linux-gnu/libldap_r-2.4.so.2 \
                    /usr/lib/x86_64-linux-gnu/liblber-2.4.so.2 \
                    /usr/lib/x86_64-linux-gnu/libsasl2.so.2 \
                    /usr/lib/x86_64-linux-gnu/libgssapi.so.3 \
                    /usr/lib/x86_64-linux-gnu/libheimntlm.so.0 \
                    /usr/lib/x86_64-linux-gnu/libkrb5.so.26 \
                    /usr/lib/x86_64-linux-gnu/libasn1.so.8 \
                    /usr/lib/x86_64-linux-gnu/libhcrypto.so.4 \
                    /usr/lib/x86_64-linux-gnu/libroken.so.18 \
                    /usr/lib/x86_64-linux-gnu/libwind.so.0 \
                    /usr/lib/x86_64-linux-gnu/libheimbase.so.1 \
                    /usr/lib/x86_64-linux-gnu/libhx509.so.5 \
                    /lib/x86_64-linux-gnu/
COPY --from=builder /usr/bin/samtools /usr/bin/
COPY --from=builder /opt/liquidator/bamliquidator \
                    /opt/liquidator/bamliquidator_bins \
                    /opt/liquidator/bamliquidator_regions \
                    ./
COPY --from=builder /opt/liquidator/bamliquidatorbatch /opt/liquidator/bamliquidatorbatch

ENV PATH="$PATH:/opt/liquidator"

RUN python3 bamliquidatorbatch/test.py

ARG GIT_COMMIT
LABEL git_commit=$GIT_COMMIT

ENTRYPOINT ["/usr/bin/python3","/opt/liquidator/bamliquidatorbatch/bamliquidator_batch.py"]
