make 
make install
COPY src/interface/phyloacc_interface.py %PREFIX%/bin/.
mkdir -p %SP_DIR%
COPY -R src/interface/phyloacc_lib %SP_DIR%/.
