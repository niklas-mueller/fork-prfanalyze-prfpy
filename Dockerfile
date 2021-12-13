FROM garikoitz/prfanalyze-base:1.0.0

RUN . /opt/conda/etc/profile.d/conda.sh \
 && conda create -n prfpy_analysis --channel intel intelpython3_full \
 && conda activate prfpy_analysis \
 && pip install --upgrade nilearn \
 && pip install --upgrade nibabel \
 && pip install --upgrade h5py \
 && pip install --upgrade wget \
 && pip install --upgrade bids \
 && pip install --upgrade sharedmem \
 && pip install --upgrade pimms \
 && pip install -U setuptools wheel \
 && git clone https://github.com/VU-Cog-Sci/prfpy.git \
 && cd prfpy \
 && python setup.py install \
 && cd .. \
 && git clone https://github.com/gallantlab/pycortex.git \
 && cd pycortex \
 && python setup.py install  \
 && cd .. \
 && mkdir -p /root/.config/pycortex/ \
 && touch /root/.config/pycortex/options.cfg 


COPY default_config.json /opt/default_config.json
COPY solve.sh /solve.sh
COPY run_prfpy.py /scripts/run_prfpy.py

RUN chmod 755 /solve.sh
# RUN chmod 755 /scripts/run_prfpy.py

ENV PRF_SOLVER prfpy